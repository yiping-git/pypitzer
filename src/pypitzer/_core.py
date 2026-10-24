# -*- coding: utf-8 -*-
# Author: Yiping Liu
# Description: Pitzer model for direct solution property calculation.
#              Accepts pre-resolved ion molalities + pre-fetched parameters.
#              No speciation, no equilibrium solving, no database access.
#              Input:  T, P, ion molalities (fully dissociated)
#              Output: activity coefficients, osmotic coefficient, water activity
# Version: 2.0

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from numpy.typing import NDArray

from . import _equations as pm


# ── Exceptions ────────────────────────────────────────────────────────────────


class PitzerError(Exception):
    pass


class MissingParameterError(PitzerError):
    def __init__(self, ion_i: "Ion", ion_j: "Ion"):
        super().__init__(
            f"No binary parameter for "
            f"{ion_i.formula}(id={ion_i.species_id}, z={ion_i.charge}) × "
            f"{ion_j.formula}(id={ion_j.species_id}, z={ion_j.charge})"
        )


class OutOfRangeError(PitzerError):
    pass


# ── Input data types ──────────────────────────────────────────────────────────


@dataclass
class Ion:
    """
    A fully resolved aqueous ion ready for Pitzer calculation.
    Constructed by the repository/facade layer from formula-keyed inputs.
    """

    species_id: int  # FK to Species.id — used as lookup key
    formula: str  # e.g. "Na+", "Cl-"  — for error messages only
    charge: int  # z — used in all Pitzer equations
    molality: float  # m [mol/kg]

    @property
    def is_cation(self) -> bool:
        return self.charge > 0

    @property
    def is_anion(self) -> bool:
        return self.charge < 0

    @property
    def is_neutral(self) -> bool:
        return self.charge == 0


@dataclass
class BinaryParam:
    """
    β⁽⁰⁾, β⁽¹⁾, β⁽²⁾, Cφ for one cation–anion pair.
    Parameters are evaluated at temperature by the repository before being
    passed to this database-independent calculation layer.
    """

    cation_id: int
    anion_id: int
    beta0: float
    beta1: float
    beta2: float = 0.0  # non-zero only for 2:2 electrolytes
    cphi: float = 0.0
    t_min: Optional[float] = None
    t_max: Optional[float] = None

    def validate_temperature(self, T_celsius: float) -> None:
        if self.t_min is not None and T_celsius < self.t_min:
            raise OutOfRangeError(
                f"T={T_celsius}°C below valid range [{self.t_min}, {self.t_max}] "
                f"for cation_id={self.cation_id}, anion_id={self.anion_id}"
            )
        if self.t_max is not None and T_celsius > self.t_max:
            raise OutOfRangeError(
                f"T={T_celsius}°C above valid range [{self.t_min}, {self.t_max}] "
                f"for cation_id={self.cation_id}, anion_id={self.anion_id}"
            )


@dataclass
class ThetaParam:
    """
    θ — like-charge interaction: cation–cation or anion–anion.
    Canonical ordering: species_i_id < species_j_id.
    """

    species_i_id: int
    species_j_id: int
    theta: float


@dataclass
class PsiParam:
    """
    ψ — triplet: two like-charge ions + one unlike-charge ion.
    Either cation–cation–anion (cca) or anion–anion–cation (aac).
    """

    psi: float
    # cca form
    cation_i_id: Optional[int] = None
    cation_j_id: Optional[int] = None
    anion_id: Optional[int] = None
    # aac form
    anion_i_id: Optional[int] = None
    anion_j_id: Optional[int] = None
    cation_id: Optional[int] = None

    @property
    def is_cca(self) -> bool:
        return self.cation_i_id is not None

    def key(self) -> tuple:
        """Canonical tuple key for indexing — like-charge pair sorted."""
        if self.is_cca:
            ci, cj = sorted([self.cation_i_id, self.cation_j_id])
            return (ci, cj, self.anion_id)
        else:
            ai, aj = sorted([self.anion_i_id, self.anion_j_id])
            return (ai, aj, self.cation_id)


@dataclass
class LambdaParam:
    """λ — neutral–ion interaction."""

    neutral_id: int
    ion_id: int
    lambda_: float


@dataclass
class ZetaParam:
    """ζ — neutral–cation–anion triplet."""

    neutral_id: int
    cation_id: int
    anion_id: int
    zeta: float


@dataclass
class PitzerInput:
    """
    Complete, self-contained input for one Pitzer property calculation.
    All parameters are pre-fetched by the repository — the model has
    zero DB awareness.

    ions            — fully dissociated aqueous ions with molalities
    binary_params   — β⁽⁰⁾, β⁽¹⁾, β⁽²⁾, Cφ for every cation–anion pair
    theta_params    — θ for like-charge pairs (optional, default 0)
    psi_params      — ψ for triplets (optional, default 0)
    lambda_params   — λ for neutral–ion pairs (optional)
    zeta_params     — ζ for neutral–cation–anion triplets (optional)
    T               — temperature [K]
    P               — pressure [MPa], None = saturation pressure
    """

    ions: list[Ion]
    binary_params: list[BinaryParam]
    theta_params: list[ThetaParam] = field(default_factory=list)
    psi_params: list[PsiParam] = field(default_factory=list)
    lambda_params: list[LambdaParam] = field(default_factory=list)
    zeta_params: list[ZetaParam] = field(default_factory=list)
    T: float = 298.15  # [K]
    P: Optional[float] = None  # [MPa]


@dataclass
class PitzerResult:
    """Output from one Pitzer calculation."""

    ionic_strength: float
    osmotic_coefficient: float
    water_activity: float
    # per-ion: keyed by species_id
    ln_gamma: dict[int, float]
    # per cation–anion pair: keyed by (cation_id, anion_id)
    mean_activity_coeff: dict[tuple[int, int], float]

    def gamma(self, species_id: int) -> float:
        """Activity coefficient γ for a single ion."""
        return float(np.exp(self.ln_gamma[species_id]))

    def log_gamma(self, species_id: int) -> float:
        """log₁₀ γ for a single ion."""
        return self.ln_gamma[species_id] / np.log(10)


# ── Main model class ──────────────────────────────────────────────────────────


class PitzerProperties:
    """
    Pure Pitzer property calculator.

    Given fully resolved ion molalities and pre-fetched parameters,
    computes activity coefficients, osmotic coefficient, and water activity.

    No speciation. No equilibrium solving. No database access.
    No string-keyed parameter lookup. Ions are identified by species_id (int).

    Usage:
        model = PitzerProperties(input)
        result = model.calculate()
    """

    # Pitzer universal constants
    ALPHA = 2.0  # kg^(1/2) mol^(-1/2)  — standard for 1:1, 2:1 etc.
    ALPHA_1 = 1.4  # kg^(1/2) mol^(-1/2)  — for 2:2 (divalent–divalent)
    ALPHA_2 = 12.0  # kg^(1/2) mol^(-1/2)  — for 2:2
    MW_WATER = 55.50844  # mol/kg

    def __init__(self, input: PitzerInput):
        self.input = input
        self.T = input.T
        self.t = input.T - 273.15  # Celsius, kept for any T-dependent calls
        self.P = input.P

        # ── Separate ions by charge ───────────────────────────────────────
        self.cations = [ion for ion in input.ions if ion.is_cation]
        self.anions = [ion for ion in input.ions if ion.is_anion]
        self.neutrals = [ion for ion in input.ions if ion.is_neutral]
        self.all_ions = input.ions

        # ── Index ions by species_id ──────────────────────────────────────
        # species_id is the canonical identity — same int that came from DB
        self.ion_index: dict[int, Ion] = {ion.species_id: ion for ion in input.ions}

        # ── Build integer position arrays for vectorised numpy ops ────────
        # position in self.all_ions list — used to build N×N matrices
        self._pos: dict[int, int] = {
            ion.species_id: i for i, ion in enumerate(self.all_ions)
        }
        self.N = len(self.all_ions)

        # charge array aligned with self.all_ions
        self.charges = np.array([ion.charge for ion in self.all_ions], dtype=np.float64)

        # molality array aligned with self.all_ions
        self.molalities = np.array(
            [ion.molality for ion in self.all_ions], dtype=np.float64
        )

        # ── Build O(1) parameter lookup indexes ───────────────────────────
        # keyed by (cation_id, anion_id) — matches BinaryParam fields
        self._binary_index: dict[tuple[int, int], BinaryParam] = {
            (p.cation_id, p.anion_id): p for p in input.binary_params
        }

        # keyed by canonical (min_id, max_id)
        self._theta_index: dict[tuple[int, int], ThetaParam] = {
            self._canon2(p.species_i_id, p.species_j_id): p for p in input.theta_params
        }

        # keyed by PsiParam.key() — like-charge pair sorted + unlike appended
        self._psi_index: dict[tuple, PsiParam] = {p.key(): p for p in input.psi_params}

        # keyed by (neutral_id, ion_id)
        self._lambda_index: dict[tuple[int, int], LambdaParam] = {
            (p.neutral_id, p.ion_id): p for p in input.lambda_params
        }

        # keyed by (neutral_id, cation_id, anion_id)
        self._zeta_index: dict[tuple[int, int, int], ZetaParam] = {
            (p.neutral_id, p.cation_id, p.anion_id): p for p in input.zeta_params
        }

        # ── Build N×N / N×N×N parameter matrices once ─────────────────────
        # These are T-independent after the repository has evaluated
        # temperature dependence into scalar beta0, beta1 etc.
        self._B: NDArray | None = None  # β matrix (ionic-strength dependent)
        self._Bphi: NDArray | None = None
        self._Bprime: NDArray | None = None
        self._C: NDArray | None = None  # Cφ/(2√|z+z-|)
        self._PHI: NDArray | None = None
        self._PHIprime: NDArray | None = None
        self._PSI: NDArray | None = None
        self._LAMBDA: NDArray | None = None
        self._ZETA: NDArray | None = None

        # ── Scalar cache (molality-dependent computed values) ─────────────
        self._cache: dict = {}

        # ── Validate binary parameter coverage and T ranges ───────────────
        self._validate()

    # ── Validation ────────────────────────────────────────────────────────────

    def _validate(self) -> None:
        """
        Every cation–anion pair MUST have a binary parameter.
        θ and ψ are optional (treated as zero when absent).
        """
        T_celsius = self.t
        for c in self.cations:
            for a in self.anions:
                param = self._binary_index.get((c.species_id, a.species_id))
                if param is None:
                    raise MissingParameterError(c, a)
                param.validate_temperature(T_celsius)

    # ── Static helpers ────────────────────────────────────────────────────────

    @staticmethod
    def _canon2(i: int, j: int) -> tuple[int, int]:
        """Canonical ordering for symmetric pairs."""
        return (i, j) if i < j else (j, i)

    def _pos_of(self, ion: Ion) -> int:
        return self._pos[ion.species_id]

    # ── Core scalar properties ────────────────────────────────────────────────

    @property
    def ionic_strength(self) -> float:
        if "I" not in self._cache:
            self._cache["I"] = float(0.5 * np.sum(self.molalities * self.charges**2))
        return self._cache["I"]

    @property
    def Z(self) -> float:
        """Z = Σ mᵢ |zᵢ|"""
        if "Z" not in self._cache:
            self._cache["Z"] = float(self.molalities @ np.abs(self.charges))
        return self._cache["Z"]

    @property
    def a_phi(self) -> float:
        """Debye–Hückel A_φ at temperature T."""
        if "a_phi" not in self._cache:
            self._cache["a_phi"] = pm.a_phi_spencer(self.T)
        return self._cache["a_phi"]

    # ── β matrices (ionic-strength dependent) ─────────────────────────────────

    def _is_divalent_pair(self, c: Ion, a: Ion) -> bool:
        return abs(c.charge) == 2 and abs(a.charge) == 2

    def _build_beta_matrices(self) -> None:
        """Build B, B_phi, B_prime, C matrices from binary parameters."""
        N = self.N
        I = self.ionic_strength
        sqrtI = np.sqrt(I)

        B = np.zeros((N, N))
        B_phi = np.zeros((N, N))
        B_prime = np.zeros((N, N))
        C = np.zeros((N, N))

        for c in self.cations:
            for a in self.anions:
                p = self._binary_index[(c.species_id, a.species_id)]
                i = self._pos_of(c)
                j = self._pos_of(a)
                b0 = p.beta0
                b1 = p.beta1
                b2 = p.beta2
                cphi = p.cphi

                if self._is_divalent_pair(c, a):
                    al1, al2 = self.ALPHA_1, self.ALPHA_2
                    g1 = pm.g_func(al1 * sqrtI)
                    g2 = pm.g_func(al2 * sqrtI)
                    gp1 = pm.g_func_prime(al1 * sqrtI)
                    gp2 = pm.g_func_prime(al2 * sqrtI)

                    b_val = b0 + b1 * g1 + b2 * g2
                    bphi_val = (
                        b0 + b1 * np.exp(-al1 * sqrtI) + b2 * np.exp(-al2 * sqrtI)
                    )
                    if I > 0:
                        bprime_val = (b1 * gp1 + b2 * gp2) / I
                    else:
                        bprime_val = 0.0
                else:
                    al = self.ALPHA
                    g = pm.g_func(al * sqrtI)
                    gp = pm.g_func_prime(al * sqrtI)

                    b_val = b0 + b1 * g
                    bphi_val = b0 + b1 * np.exp(-al * sqrtI)
                    if I > 0:
                        bprime_val = b1 * gp / I
                    else:
                        bprime_val = 0.0

                c_val = cphi / (2.0 * np.sqrt(abs(c.charge * a.charge)))

                # symmetric
                B[i, j] = B[j, i] = b_val
                B_phi[i, j] = B_phi[j, i] = bphi_val
                B_prime[i, j] = B_prime[j, i] = bprime_val
                C[i, j] = C[j, i] = c_val

        self._B = B
        self._Bphi = B_phi
        self._Bprime = B_prime
        self._C = C

    @property
    def B(self) -> NDArray:
        if self._B is None:
            self._build_beta_matrices()
        return self._B

    @property
    def Bphi(self) -> NDArray:
        if self._Bphi is None:
            self._build_beta_matrices()
        return self._Bphi

    @property
    def Bprime(self) -> NDArray:
        if self._Bprime is None:
            self._build_beta_matrices()
        return self._Bprime

    @property
    def C_matrix(self) -> NDArray:
        if self._C is None:
            self._build_beta_matrices()
        return self._C

    # ── Φ matrices (like-charge, I-dependent via e_theta) ─────────────────────

    def _build_phi_matrices(self) -> None:
        N = self.N
        I = self.ionic_strength
        Af = self.a_phi
        PHI = np.zeros((N, N))
        PHIprime = np.zeros((N, N))

        # gather all like-charge pairs that have a theta param
        for p in self.input.theta_params:
            ion_i = self.ion_index[p.species_i_id]
            ion_j = self.ion_index[p.species_j_id]
            i = self._pos_of(ion_i)
            j = self._pos_of(ion_j)
            z1, z2 = ion_i.charge, ion_j.charge

            e_theta, e_theta_prime = pm.get_e_theta(z1, z2, Af, I)
            phi_val = p.theta + e_theta
            phi_prime_val = e_theta_prime

            PHI[i, j] = PHI[j, i] = phi_val
            PHIprime[i, j] = PHIprime[j, i] = phi_prime_val

        self._PHI = PHI
        self._PHIprime = PHIprime

    @property
    def PHI(self) -> NDArray:
        if self._PHI is None:
            self._build_phi_matrices()
        return self._PHI

    @property
    def PHIprime(self) -> NDArray:
        if self._PHIprime is None:
            self._build_phi_matrices()
        return self._PHIprime

    @property
    def PHI_phi(self) -> NDArray:
        """Φφ = Φ + I·Φ' — used in osmotic coefficient."""
        return self.PHI + self.ionic_strength * self.PHIprime

    # ── Ψ matrix (triplets, T-only dependent) ─────────────────────────────────

    @property
    def PSI(self) -> NDArray:
        if self._PSI is None:
            N = self.N
            PSI = np.zeros((N, N, N))
            for p in self.input.psi_params:
                if p.is_cca:
                    ci = self._pos[p.cation_i_id]
                    cj = self._pos[p.cation_j_id]
                    a = self._pos[p.anion_id]
                    PSI[ci, cj, a] = PSI[cj, ci, a] = p.psi
                else:
                    ai = self._pos[p.anion_i_id]
                    aj = self._pos[p.anion_j_id]
                    c = self._pos[p.cation_id]
                    PSI[ai, aj, c] = PSI[aj, ai, c] = p.psi
            self._PSI = PSI
        return self._PSI

    # ── Λ / ζ matrices (neutral interactions) ────────────────────────────────

    @property
    def LAMBDA(self) -> NDArray:
        if self._LAMBDA is None:
            N = self.N
            L = np.zeros((N, N))
            for p in self.input.lambda_params:
                n = self._pos[p.neutral_id]
                i = self._pos[p.ion_id]
                L[n, i] = L[i, n] = p.lambda_
            self._LAMBDA = L
        return self._LAMBDA

    @property
    def ZETA(self) -> NDArray:
        if self._ZETA is None:
            N = self.N
            Z = np.zeros((N, N, N))
            for p in self.input.zeta_params:
                n = self._pos[p.neutral_id]
                c = self._pos[p.cation_id]
                a = self._pos[p.anion_id]
                Z[n, c, a] = Z[n, a, c] = p.zeta
            self._ZETA = Z
        return self._ZETA

    # ── F function ────────────────────────────────────────────────────────────

    @property
    def F(self) -> float:
        """
        F = f_γ + Σ_ca m_c m_a B'_ca + Σ_cc m_c m_c Φ'_cc + Σ_aa m_a m_a Φ'_aa
        """
        if "F" in self._cache:
            return self._cache["F"]

        m = self.molalities
        I = self.ionic_strength
        Af = self.a_phi

        f_gamma = pm.get_f_gamma(Af, I)

        c_idx = [self._pos_of(ion) for ion in self.cations]
        a_idx = [self._pos_of(ion) for ion in self.anions]

        val = f_gamma

        # ca term
        if c_idx and a_idx:
            Bp = self.Bprime[np.ix_(c_idx, a_idx)]
            val += np.einsum("i,j,ij->", m[c_idx], m[a_idx], Bp)

        # cc term
        cc_pairs = self._like_charge_pairs(self.cations)
        for i, j in cc_pairs:
            val += m[i] * m[j] * self.PHIprime[i, j]

        # aa term
        aa_pairs = self._like_charge_pairs(self.anions)
        for i, j in aa_pairs:
            val += m[i] * m[j] * self.PHIprime[i, j]

        self._cache["F"] = float(val)
        return self._cache["F"]

    # ── Activity coefficients ─────────────────────────────────────────────────

    def _ln_gamma_ion(self, ion: Ion) -> float:
        """
        ln γ for a single ion (cation or anion).
        Full Pitzer expression including β, C, Φ, Ψ, Λ, ζ terms.
        """
        m = self.molalities
        z = ion.charge
        pos = self._pos_of(ion)

        c_idx = np.array([self._pos_of(c) for c in self.cations])
        a_idx = np.array([self._pos_of(a) for a in self.anions])
        n_idx = np.array([self._pos_of(n) for n in self.neutrals])

        m_c = m[c_idx]
        m_a = m[a_idx]

        # term 0: z² F
        t0 = z**2 * self.F

        # term 1: Σ_a m_a (2B_Ma + Z C_Ma)  for cation M
        #         Σ_c m_c (2B_Xc + Z C_Xc)  for anion X
        if ion.is_cation:
            B_row = self.B[pos, a_idx]
            C_row = self.C_matrix[pos, a_idx]
            t1 = float(np.dot(m_a, 2.0 * B_row + self.Z * C_row))
        else:
            B_row = self.B[pos, c_idx]
            C_row = self.C_matrix[pos, c_idx]
            t1 = float(np.dot(m_c, 2.0 * B_row + self.Z * C_row))

        # term 2: like-charge Φ + Ψ sums
        t2 = 0.0
        if ion.is_cation:
            for other in self.cations:
                if other.species_id == ion.species_id:
                    continue
                oi = self._pos_of(other)
                phi_val = self.PHI[pos, oi]
                psi_sum = sum(
                    m[self._pos_of(a)] * self.PSI[pos, oi, self._pos_of(a)]
                    for a in self.anions
                )
                t2 += m[oi] * (2.0 * phi_val + psi_sum)
        else:
            for other in self.anions:
                if other.species_id == ion.species_id:
                    continue
                oi = self._pos_of(other)
                phi_val = self.PHI[pos, oi]
                psi_sum = sum(
                    m[self._pos_of(c)] * self.PSI[pos, oi, self._pos_of(c)]
                    for c in self.cations
                )
                t2 += m[oi] * (2.0 * phi_val + psi_sum)

        # term 3: Ψ from unlike-charge pairs
        # for M cation: Σ_aa m_a1 m_a2 Ψ_a1a2M
        # for X anion:  Σ_cc m_c1 m_c2 Ψ_c1c2X
        t3 = 0.0
        if ion.is_cation:
            aa_pairs = self._like_charge_pairs(self.anions)
            for ai, aj in aa_pairs:
                t3 += m[ai] * m[aj] * self.PSI[ai, aj, pos]
        else:
            cc_pairs = self._like_charge_pairs(self.cations)
            for ci, cj in cc_pairs:
                t3 += m[ci] * m[cj] * self.PSI[ci, cj, pos]

        # term 4: |z| Σ_ca m_c m_a C_ca
        t4 = 0.0
        if c_idx.size > 0 and a_idx.size > 0:
            C_sub = self.C_matrix[np.ix_(c_idx, a_idx)]
            t4 = abs(z) * float(np.einsum("i,j,ij->", m_c, m_a, C_sub))

        # term 5: neutral λ terms  2 Σ_n m_n λ_nM
        t5 = 0.0
        if n_idx.size > 0:
            for n in self.neutrals:
                ni = self._pos_of(n)
                lam = self.LAMBDA[ni, pos]
                t5 += 2.0 * m[ni] * lam

        # term 6: neutral ζ terms
        t6 = 0.0
        if n_idx.size > 0:
            for n in self.neutrals:
                ni = self._pos_of(n)
                if ion.is_cation:
                    for a in self.anions:
                        ai = self._pos_of(a)
                        t6 += m[ni] * m[ai] * self.ZETA[ni, pos, ai]
                else:
                    for c in self.cations:
                        ci = self._pos_of(c)
                        t6 += m[ni] * m[ci] * self.ZETA[ni, ci, pos]

        return t0 + t1 + t2 + t3 + t4 + t5 + t6

    def _ln_gamma_neutral(self, ion: Ion) -> float:
        """ln γ for a neutral species (e.g. CO₂)."""
        m = self.molalities
        pos = self._pos_of(ion)

        t0 = sum(
            2.0 * m[self._pos_of(c)] * self.LAMBDA[pos, self._pos_of(c)]
            for c in self.cations
        )
        t1 = sum(
            2.0 * m[self._pos_of(a)] * self.LAMBDA[pos, self._pos_of(a)]
            for a in self.anions
        )
        t2 = sum(
            m[self._pos_of(c)]
            * m[self._pos_of(a)]
            * self.ZETA[pos, self._pos_of(c), self._pos_of(a)]
            for c in self.cations
            for a in self.anions
        )
        return t0 + t1 + t2

    # ── Osmotic coefficient ───────────────────────────────────────────────────

    def _compute_osmotic_coefficient(self) -> float:
        """
        φ = 1 + (2/Σm) × [DH term + B_phi term + C term + Φ_phi term + Ψ term + Λ term + ζ term]
        """
        m = self.molalities
        I = self.ionic_strength
        Af = self.a_phi
        m_sum = float(np.sum(m))
        if m_sum == 0:
            return 1.0

        c_idx = np.array([self._pos_of(c) for c in self.cations])
        a_idx = np.array([self._pos_of(a) for a in self.anions])
        n_idx = np.array([self._pos_of(n) for n in self.neutrals])

        m_c = m[c_idx]
        m_a = m[a_idx]

        # DH term
        dh = -(Af * I**1.5) / (1.0 + 1.2 * np.sqrt(I))

        # B_phi + C term — ca pairs
        bphi_c = 0.0
        if c_idx.size > 0 and a_idx.size > 0:
            Bp = self.Bphi[np.ix_(c_idx, a_idx)]
            C = self.C_matrix[np.ix_(c_idx, a_idx)]
            bphi_c = float(np.einsum("i,j,ij->", m_c, m_a, Bp + self.Z * C))

        # Φ_phi + Ψ — cc pairs
        cc_term = 0.0
        cc_pairs = self._like_charge_pairs(self.cations)
        for ci, cj in cc_pairs:
            psi_sum = sum(
                m[self._pos_of(a)] * self.PSI[ci, cj, self._pos_of(a)]
                for a in self.anions
            )
            cc_term += m[ci] * m[cj] * (self.PHI_phi[ci, cj] + psi_sum)

        # Φ_phi + Ψ — aa pairs
        aa_term = 0.0
        aa_pairs = self._like_charge_pairs(self.anions)
        for ai, aj in aa_pairs:
            psi_sum = sum(
                m[self._pos_of(c)] * self.PSI[ai, aj, self._pos_of(c)]
                for c in self.cations
            )
            aa_term += m[ai] * m[aj] * (self.PHI_phi[ai, aj] + psi_sum)

        # Λ — neutral terms
        lambda_term = 0.0
        if n_idx.size > 0:
            for n in self.neutrals:
                ni = self._pos_of(n)
                for c in self.cations:
                    ci = self._pos_of(c)
                    lambda_term += m[ni] * m[ci] * self.LAMBDA[ni, ci]
                for a in self.anions:
                    ai = self._pos_of(a)
                    lambda_term += m[ni] * m[ai] * self.LAMBDA[ni, ai]

        # ζ — nca triplets
        zeta_term = 0.0
        if n_idx.size > 0:
            for n in self.neutrals:
                ni = self._pos_of(n)
                for c in self.cations:
                    ci = self._pos_of(c)
                    for a in self.anions:
                        ai = self._pos_of(a)
                        zeta_term += m[ni] * m[ci] * m[ai] * self.ZETA[ni, ci, ai]

        inner = dh + bphi_c + cc_term + aa_term + lambda_term + zeta_term
        return 1.0 + (2.0 / m_sum) * inner

    # ── Water activity ────────────────────────────────────────────────────────

    def _compute_water_activity(self, phi: float) -> float:
        """ln(a_w) = -φ · Σm / 55.508"""
        m_sum = float(np.sum(self.molalities))
        return -phi * m_sum / self.MW_WATER

    # ── Mean activity coefficient ─────────────────────────────────────────────

    def _mean_activity_coeff(
        self, cation: Ion, anion: Ion, ln_gamma: dict[int, float]
    ) -> float:
        """
        γ± = (γ+^ν+ · γ-^ν-)^(1/ν)
        For a simple 1:1 salt ν+=ν-=1 → γ± = √(γ+ · γ-)
        """
        nu_c = abs(anion.charge)  # stoichiometry from charge balance
        nu_a = abs(cation.charge)
        nu = nu_c + nu_a
        ln_gpm = (
            nu_c * ln_gamma[cation.species_id] + nu_a * ln_gamma[anion.species_id]
        ) / nu
        return float(np.exp(ln_gpm))

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _like_charge_pairs(self, ions: list[Ion]) -> list[tuple[int, int]]:
        """
        All unique pairs from a list of ions (same charge sign).
        Returns position-index tuples for direct numpy indexing.
        """
        pairs = []
        for i in range(len(ions)):
            for j in range(i + 1, len(ions)):
                pairs.append((self._pos_of(ions[i]), self._pos_of(ions[j])))
        return pairs

    # ── Main entry point ──────────────────────────────────────────────────────

    def calculate(self) -> PitzerResult:
        """
        Run the full property calculation.
        Returns PitzerResult with all solution properties.
        """
        # activity coefficients — one per ion
        ln_gamma: dict[int, float] = {}
        for ion in self.all_ions:
            if ion.is_neutral:
                ln_gamma[ion.species_id] = self._ln_gamma_neutral(ion)
            else:
                ln_gamma[ion.species_id] = self._ln_gamma_ion(ion)

        # osmotic coefficient
        phi = self._compute_osmotic_coefficient()

        # water activity
        ln_aw = self._compute_water_activity(phi)

        # mean activity coefficients for every cation–anion pair
        mean_gamma: dict[tuple[int, int], float] = {}
        for c in self.cations:
            for a in self.anions:
                mean_gamma[(c.species_id, a.species_id)] = self._mean_activity_coeff(
                    c, a, ln_gamma
                )

        return PitzerResult(
            ionic_strength=self.ionic_strength,
            osmotic_coefficient=phi,
            water_activity=float(np.exp(ln_aw)),
            ln_gamma=ln_gamma,
            mean_activity_coeff=mean_gamma,
        )

