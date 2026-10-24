"""Read-only access to the bundled, versioned thermodynamic database."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from itertools import combinations, product
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


def canonical_reaction(equation: str) -> str:
    try:
        lhs, rhs = equation.split("=", 1)
    except ValueError as error:
        raise ValueError("Reaction must contain '='") from error
    split_side = lambda side: sorted(part.strip() for part in re.split(r"\s+\+\s+", side.strip()) if part.strip())
    return f"{' + '.join(split_side(lhs))} = {' + '.join(split_side(rhs))}"


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

    def reaction(self, equation: str, temperature_k: float) -> Reaction:
        canonical = canonical_reaction(equation)
        reaction_id = hashlib.sha256(canonical.encode()).hexdigest()
        with self._connect() as connection:
            row = connection.execute(
                f"SELECT * FROM reaction WHERE reaction_id = ? AND {self._temperature_clause()} ORDER BY T_min LIMIT 1",
                (reaction_id, temperature_k, temperature_k),
            ).fetchone()
        if row is None:
            raise ReactionNotFoundError(f"No reaction data for {equation!r} at {temperature_k:g} K")
        stored_equation = str(row["reaction"])
        return Reaction(
            equation=stored_equation,
            ln_k=_evaluate(row["analytic"], temperature_k, row["eq_num"]),
            stoichiometry=_reaction_stoichiometry(stored_equation),
            temperature_min_k=float(row["T_min"]),
            temperature_max_k=float(row["T_max"]),
            reference=row["ref"],
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
