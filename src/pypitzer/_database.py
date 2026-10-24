"""Read-only access to the bundled, versioned thermodynamic database."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
import json
from pathlib import Path
import re
import sqlite3
from types import MappingProxyType
from typing import Mapping

from . import _equations
from ._core import BinaryParam, Ion, LambdaParam, PitzerInput, PsiParam, ThetaParam, ZetaParam

DATABASE_VERSION = "2026.1"
DATABASE_SHA256 = "13a5c8ac3d6323218898ca85758f3ff67ba3235cd90148c492e4abf7aae2eadf"
DEFAULT_DATABASE = Path(__file__).with_name("data") / "parameters-v2026.1.sqlite"


class DatabaseError(RuntimeError):
    """Base class for runtime database errors."""


class ParameterNotFoundError(DatabaseError):
    """Raised when required parameter data are absent."""


class ReactionNotFoundError(DatabaseError):
    """Raised when a reaction is absent or invalid at the requested temperature."""


@dataclass(frozen=True)
class Reaction:
    equation: str
    ln_k: float
    stoichiometry: Mapping[str, float]
    temperature_min_k: float
    temperature_max_k: float
    reference: str | None = None


@dataclass(frozen=True)
class PhaseRecord:
    phase_id: str
    formula: str
    temperature_min_k: float
    temperature_max_k: float


_EXCLUDED_PHASE_REACTIONS = {
    # This record has solid species on both sides and is not a dissolution
    # equilibrium suitable for the fluid-inclusion solver.
    "b96c269e3255eeaa1ce21ce25ac454a029774c2ae8d2a88ebe80db69b7e4a0b7",
}

_PHASE_ALIASES = {
    "ice": "H2O",
    "halite": "NaCl",
    "sylvite": "KCl",
    "hydrohalite": "NaCl-2H2O",
}

def charge_of(formula: str) -> int:
    """Return the charge from geochemical labels such as ``Mg+2`` or ``Cl-``."""
    species = re.sub(r"\((aq|s|l|g)\)$", "", formula.strip(), flags=re.IGNORECASE)
    match = re.search(r"([+-])(\d*)$", species)
    if match is None:
        return 0
    magnitude = int(match.group(2) or "1")
    return magnitude if match.group(1) == "+" else -magnitude


def strip_state(formula: str) -> str:
    return re.sub(r"\((aq|s|l|g)\)$", "", formula.strip(), flags=re.IGNORECASE)


def canonical_group(species: tuple[str, ...] | list[str]) -> str:
    return ",".join(sorted(strip_state(item) for item in species))


def _coefficients(text: str | None) -> list[float]:
    return [] if not text else [float(value) for value in text.split(",")]


def _evaluate(text: str | None, temperature_k: float, equation: int | None) -> float:
    values = _coefficients(text)
    if not values:
        return 0.0
    return _equations.evaluate_coefficients(values, temperature_k, int(equation or 0))


def _parse_term(term: str) -> tuple[str, float]:
    match = re.match(r"^\s*(?:(\d+(?:\.\d+)?)\s*)?(.+?)\s*$", term)
    if match is None:
        raise ValueError(f"Invalid reaction term {term!r}")
    return match.group(2), float(match.group(1) or 1.0)


def _reaction_stoichiometry(equation: str) -> Mapping[str, float]:
    lhs, rhs = equation.split("=", 1)
    result: dict[str, float] = {}
    for sign, side in ((-1.0, lhs), (1.0, rhs)):
        for token in re.split(r"\s+\+\s+", side.strip()):
            formula, coefficient = _parse_term(token)
            result[formula.strip()] = result.get(formula.strip(), 0.0) + sign * coefficient
    return MappingProxyType(result)


class ParameterRepository:
    """A small SQLite repository opened permanently in read-only mode."""

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path is not None else DEFAULT_DATABASE
        if not self.path.is_file():
            raise DatabaseError(f"Thermodynamic database does not exist: {self.path}")

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(f"file:{self.path.resolve()}?mode=ro", uri=True)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        return connection

    @staticmethod
    def _temperature_clause() -> str:
        # Preserve the database's published interval convention. At a shared
        # boundary, the lower-temperature record owns the endpoint.
        return "T_min < ? AND T_max >= ?"

    def _binary_row(self, group: tuple[str, str], temperature_k: float, *, required: bool) -> sqlite3.Row | None:
        key = canonical_group(group)
        with self._connect() as connection:
            row = connection.execute(
                f"SELECT * FROM binary_params WHERE pair = ? AND {self._temperature_clause()} ORDER BY T_min LIMIT 1",
                (key, temperature_k, temperature_k),
            ).fetchone()
        if row is None and required:
            raise ParameterNotFoundError(f"No binary parameter for {key} at {temperature_k:g} K")
        return row

    def _ternary_row(self, group: tuple[str, str, str], temperature_k: float) -> sqlite3.Row | None:
        key = canonical_group(group)
        if "H+" in group and "OH-" in group:
            key = "H+,OH-,X"
        with self._connect() as connection:
            return connection.execute(
                f"SELECT * FROM ternary_params WHERE triplet = ? AND {self._temperature_clause()} ORDER BY T_min LIMIT 1",
                (key, temperature_k, temperature_k),
            ).fetchone()

    @staticmethod
    def _reaction_from_row(row: sqlite3.Row, temperature_k: float) -> Reaction:
        stored_equation = str(row["reaction"])
        return Reaction(
            equation=stored_equation,
            ln_k=_evaluate(row["analytic"], temperature_k, row["eq_num"]),
            stoichiometry=_reaction_stoichiometry(stored_equation),
            temperature_min_k=float(row["T_min"]),
            temperature_max_k=float(row["T_max"]),
            reference=row["ref"],
        )

    @staticmethod
    def _phase_identity(row: sqlite3.Row) -> tuple[str, str]:
        formula = strip_state(str(row["reaction"]).split("=", 1)[0].strip())
        polymorph_data = json.loads(row["polymorph"]) if row["polymorph"] else None
        if isinstance(polymorph_data, dict) and polymorph_data:
            phase_id = str(next(iter(polymorph_data.values()))).lower()
        else:
            phase_id = formula
        return phase_id, formula

    def phases(self) -> tuple[PhaseRecord, ...]:
        """Return the curated phases available to fluid-inclusion calculations."""
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT * FROM reaction WHERE solubility = 1 ORDER BY reaction, T_min"
            ).fetchall()
        grouped: dict[str, dict[str, object]] = {}
        for row in rows:
            if row["reaction_id"] in _EXCLUDED_PHASE_REACTIONS:
                continue
            phase_id, formula = self._phase_identity(row)
            entry = grouped.setdefault(
                phase_id,
                {
                    "formula": formula,
                    "minimum": float(row["T_min"]),
                    "maximum": float(row["T_max"]),
                },
            )
            entry["minimum"] = min(float(entry["minimum"]), float(row["T_min"]))
            entry["maximum"] = max(float(entry["maximum"]), float(row["T_max"]))
        return tuple(
            PhaseRecord(
                phase_id=phase_id,
                formula=str(data["formula"]),
                temperature_min_k=float(data["minimum"]),
                temperature_max_k=float(data["maximum"]),
            )
            for phase_id, data in sorted(grouped.items(), key=lambda item: item[0].lower())
        )

    def phase_reaction(self, phase: str, temperature_k: float) -> Reaction:
        """Resolve a stable phase ID, formula, or supported alias to a reaction."""
        requested = phase.strip()
        canonical = _PHASE_ALIASES.get(requested.lower(), requested)
        known = {record.phase_id.lower(): record.phase_id for record in self.phases()}
        resolved_id = known.get(canonical.lower())
        if resolved_id is None:
            available = ", ".join(record.phase_id for record in self.phases())
            raise ReactionNotFoundError(
                f"Unknown last-melting phase {phase!r}. Supported phases: {available}"
            )
        with self._connect() as connection:
            rows = connection.execute(
                f"SELECT * FROM reaction WHERE solubility = 1 AND {self._temperature_clause()} ORDER BY T_min",
                (temperature_k, temperature_k),
            ).fetchall()
        for row in rows:
            if row["reaction_id"] in _EXCLUDED_PHASE_REACTIONS:
                continue
            phase_id, _ = self._phase_identity(row)
            if phase_id == resolved_id:
                return self._reaction_from_row(row, temperature_k)
        raise ReactionNotFoundError(
            f"Last-melting phase {resolved_id!r} has no data at {temperature_k:g} K"
        )

    def pitzer_input(self, molalities: Mapping[str, float], temperature_k: float) -> tuple[PitzerInput, Mapping[str, int]]:
        """Resolve formula-keyed molalities into a complete typed calculator input."""
        if temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive")
        clean = {strip_state(name): float(value) for name, value in molalities.items()}
        if not clean or any(value < 0.0 for value in clean.values()):
            raise ValueError("molalities must be a non-empty mapping of non-negative values")
        formulas = sorted(clean)
        ids = MappingProxyType({formula: index + 1 for index, formula in enumerate(formulas)})
        ions = [Ion(ids[name], name, charge_of(name), clean[name]) for name in formulas]
        cations = [name for name in formulas if charge_of(name) > 0]
        anions = [name for name in formulas if charge_of(name) < 0]
        neutrals = [name for name in formulas if charge_of(name) == 0]

        binary: list[BinaryParam] = []
        for cation, anion in product(cations, anions):
            row = self._binary_row((cation, anion), temperature_k, required=True)
            assert row is not None
            binary.append(BinaryParam(ids[cation], ids[anion], _evaluate(row["b0"], temperature_k, row["eq_num"]), _evaluate(row["b1"], temperature_k, row["eq_num"]), _evaluate(row["b2"], temperature_k, row["eq_num"]), _evaluate(row["c_phi"], temperature_k, row["eq_num"]), float(row["T_min"] - 273.15), float(row["T_max"] - 273.15)))

        theta: list[ThetaParam] = []
        for pair in list(combinations(cations, 2)) + list(combinations(anions, 2)):
            row = self._binary_row(pair, temperature_k, required=False)
            value = 0.0 if row is None else _evaluate(row["theta"], temperature_k, row["eq_num"])
            theta.append(ThetaParam(ids[pair[0]], ids[pair[1]], value))

        psi: list[PsiParam] = []
        for c1, c2 in combinations(cations, 2):
            for anion in anions:
                row = self._ternary_row((c1, c2, anion), temperature_k)
                if row is not None:
                    psi.append(PsiParam(_evaluate(row["psi"], temperature_k, row["eq_num"]), cation_i_id=ids[c1], cation_j_id=ids[c2], anion_id=ids[anion]))
        for a1, a2 in combinations(anions, 2):
            for cation in cations:
                row = self._ternary_row((a1, a2, cation), temperature_k)
                if row is not None:
                    psi.append(PsiParam(_evaluate(row["psi"], temperature_k, row["eq_num"]), anion_i_id=ids[a1], anion_j_id=ids[a2], cation_id=ids[cation]))

        lambdas: list[LambdaParam] = []
        for neutral, ion in product(neutrals, cations + anions):
            row = self._binary_row((neutral, ion), temperature_k, required=False)
            if row is not None:
                lambdas.append(LambdaParam(ids[neutral], ids[ion], _evaluate(row["lambda"], temperature_k, row["eq_num"])))
        zetas: list[ZetaParam] = []
        for neutral, cation, anion in product(neutrals, cations, anions):
            row = self._ternary_row((neutral, cation, anion), temperature_k)
            if row is not None:
                zetas.append(ZetaParam(ids[neutral], ids[cation], ids[anion], _evaluate(row["zeta"], temperature_k, row["eq_num"])))

        return PitzerInput(ions, binary, theta, psi, lambdas, zetas, temperature_k), ids
