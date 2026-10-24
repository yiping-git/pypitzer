"""Fluid-inclusion quantification from internal-standard molar ratios."""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.optimize import least_squares

from ._database import ParameterRepository, Reaction, charge_of, strip_state
from ._properties import PropertyCalculator, SolutionProperties

# The published implementation and the curated database use 273.16 for the
# Celsius-to-Kelvin conversion. Retaining that convention is required for
# exact reproduction of the published fluid-inclusion examples.
PUBLISHED_CELSIUS_OFFSET = 273.16

DEFAULT_ELEMENT_SPECIES = MappingProxyType(
    {
        "Ba": "Ba+2",
        "Ca": "Ca+2",
        "Cs": "Cs+",
        "Fe": "Fe+2",
        "K": "K+",
        "Li": "Li+",
        "Mg": "Mg+2",
        "Mn": "Mn+2",
        "Na": "Na+",
        "Pb": "Pb+2",
        "Rb": "Rb+",
        "Sr": "Sr+2",
        "Zn": "Zn+2",
    }
)


@dataclass(frozen=True)
class FluidInclusionProblem:
    """Inputs for a two-constraint internal-standard calculation.

    ``ratios_to_internal_standard`` are aqueous-species/internal-standard
    molar ratios. The standard itself is omitted; its ratio is always one.
    """

    ratios_to_internal_standard: Mapping[str, float]
    last_melting_phase: str
    temperature_c: float
    internal_standard: str = "Na+"
    balancing_ion: str = "Cl-"
    initial_molalities: tuple[float, float] = (1.0, 1.0)

    def __post_init__(self) -> None:
        if not self.last_melting_phase.strip():
            raise ValueError("last_melting_phase must not be empty")
        if charge_of(self.internal_standard) == 0 or charge_of(self.balancing_ion) == 0:
            raise ValueError("internal_standard and balancing_ion must be charged aqueous species")
        if self.internal_standard == self.balancing_ion:
            raise ValueError("balancing_ion must differ from internal_standard")
        if len(self.initial_molalities) != 2 or any(value <= 0.0 for value in self.initial_molalities):
            raise ValueError("initial_molalities must contain two positive values")
        normalized = {
            strip_state(name): float(value)
            for name, value in self.ratios_to_internal_standard.items()
        }
        if len(normalized) != len(self.ratios_to_internal_standard):
            raise ValueError(
                "ratios_to_internal_standard contains duplicate species after removing state labels"
            )
        if any(value < 0.0 for value in normalized.values()):
            raise ValueError("ratios_to_internal_standard must be non-negative")
        supplied_standard = normalized.get(strip_state(self.internal_standard), 1.0)
        if not math.isclose(float(supplied_standard), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("the internal-standard ratio must be 1")
        if strip_state(self.balancing_ion) in normalized:
            raise ValueError(
                "balancing_ion is solved and must not appear in ratios_to_internal_standard"
            )

    @classmethod
    def from_element_ratios(
        cls,
        element_ratios_to_internal_standard: Mapping[str, float],
        last_melting_phase: str,
        temperature_c: float,
        *,
        internal_standard_element: str = "Na",
        species_by_element: Mapping[str, str] | None = None,
        balancing_ion: str = "Cl-",
        initial_molalities: tuple[float, float] = (1.0, 1.0),
    ) -> "FluidInclusionProblem":
        """Convert measured element/standard molar ratios to species constraints.

        The defaults cover database-supported aqueous cations. Pass
        ``species_by_element`` to override an oxidation state or aqueous form.
        """
        mapping = dict(DEFAULT_ELEMENT_SPECIES)
        if species_by_element:
            mapping.update(species_by_element)
        unknown = sorted(set(element_ratios_to_internal_standard) - mapping.keys())
        if internal_standard_element not in mapping:
            unknown.append(internal_standard_element)
        if unknown:
            raise ValueError(f"No aqueous-species mapping for element(s): {', '.join(sorted(set(unknown)))}")
        ratios = {
            mapping[element]: float(value)
            for element, value in element_ratios_to_internal_standard.items()
        }
        return cls(
            ratios_to_internal_standard=ratios,
            last_melting_phase=last_melting_phase,
            temperature_c=temperature_c,
            internal_standard=mapping[internal_standard_element],
            balancing_ion=balancing_ion,
            initial_molalities=initial_molalities,
        )


@dataclass(frozen=True)
class SolverDiagnostics:
    converged: bool
    message: str
    iterations: int
    function_evaluations: int
    charge_balance_residual: float
    equilibrium_residual: float
    residual_norm: float


@dataclass(frozen=True)
class FluidInclusionResult:
    molalities: Mapping[str, float]
    last_melting_phase: str
    temperature_c: float
    diagnostics: SolverDiagnostics

    def require_converged(self) -> "FluidInclusionResult":
        if not self.diagnostics.converged:
            raise ConvergenceError(self.diagnostics.message, self)
        return self


class ConvergenceError(RuntimeError):
    def __init__(self, message: str, result: FluidInclusionResult | None = None):
        super().__init__(message)
        self.result = result


@dataclass(frozen=True)
class LastMeltingPhase:
    phase_id: str
    formula: str
    temperature_min_c: float
    temperature_max_c: float


def supported_last_melting_phases() -> tuple[LastMeltingPhase, ...]:
    """Return phase choices suitable for a PySILLS dropdown."""
    return tuple(
        LastMeltingPhase(
            phase_id=record.phase_id,
            formula=record.formula,
            temperature_min_c=record.temperature_min_k - PUBLISHED_CELSIUS_OFFSET,
            temperature_max_c=record.temperature_max_k - PUBLISHED_CELSIUS_OFFSET,
        )
        for record in ParameterRepository().phases()
    )


class FluidInclusionSolver:
    def __init__(self) -> None:
        self._repository = ParameterRepository()
        self._calculator = PropertyCalculator(self._repository)

    @staticmethod
    def _molalities(problem: FluidInclusionProblem, values: np.ndarray) -> dict[str, float]:
        standard_m, balancing_m = np.exp(values)
        result = {
            strip_state(name): float(ratio) * float(standard_m)
            for name, ratio in problem.ratios_to_internal_standard.items()
            if strip_state(name) != strip_state(problem.internal_standard)
        }
        result[strip_state(problem.internal_standard)] = float(standard_m)
        result[strip_state(problem.balancing_ion)] = float(balancing_m)
        return result

    @staticmethod
    def _charge_balance(molalities: Mapping[str, float]) -> float:
        return float(sum(charge_of(name) * value for name, value in molalities.items()))

    @staticmethod
    def _state(formula: str) -> str | None:
        match = re.search(r"\((aq|s|l|g)\)$", formula, flags=re.IGNORECASE)
        return None if match is None else match.group(1).lower()

    def _equilibrium_residual(self, reaction: Reaction, properties: SolutionProperties) -> float:
        ln_q = 0.0
        for labelled_formula, coefficient in reaction.stoichiometry.items():
            state = self._state(labelled_formula)
            formula = strip_state(labelled_formula)
            if state in {"s", "g"}:
                continue
            if formula == "H2O" and state == "l":
                ln_q += coefficient * math.log(properties.water_activity)
                continue
            if formula not in properties.molalities:
                raise ValueError(
                    f"Reaction species {formula!r} is not constrained by "
                    "ratios_to_internal_standard"
                )
            molality = properties.molalities[formula]
            if molality <= 0.0:
                raise ValueError(f"Reaction species {formula!r} has zero activity")
            ln_q += coefficient * (math.log(molality) + properties.ln_activity_coefficients[formula])
        return float(reaction.ln_k - ln_q)

    def solve(self, problem: FluidInclusionProblem, *, raise_on_failure: bool = True, tolerance: float = 1e-9, max_evaluations: int = 300) -> FluidInclusionResult:
        temperature_k = float(problem.temperature_c) + PUBLISHED_CELSIUS_OFFSET
        reaction = self._repository.phase_reaction(
            problem.last_melting_phase, temperature_k
        )

        def residual(log_values: np.ndarray) -> np.ndarray:
            molalities = self._molalities(problem, log_values)
            properties = self._calculator.calculate(molalities, temperature_k)
            return np.array([self._charge_balance(molalities), self._equilibrium_residual(reaction, properties)])

        optimization = least_squares(
            residual,
            np.log(np.asarray(problem.initial_molalities, dtype=float)),
            xtol=tolerance,
            ftol=tolerance,
            gtol=tolerance,
            max_nfev=max_evaluations,
        )
        molalities = self._molalities(problem, optimization.x)
        properties = self._calculator.calculate(molalities, temperature_k)
        final_residual = residual(optimization.x)
        converged = bool(optimization.success and np.max(np.abs(final_residual)) <= max(tolerance * 10.0, 1e-8))
        diagnostics = SolverDiagnostics(
            converged=converged,
            message=str(optimization.message),
            iterations=int(optimization.njev or 0),
            function_evaluations=int(optimization.nfev),
            charge_balance_residual=float(final_residual[0]),
            equilibrium_residual=float(final_residual[1]),
            residual_norm=float(np.linalg.norm(final_residual)),
        )
        result = FluidInclusionResult(
            MappingProxyType(molalities),
            problem.last_melting_phase,
            float(problem.temperature_c),
            diagnostics,
        )
        if raise_on_failure and not converged:
            raise ConvergenceError(f"Fluid-inclusion solver did not converge: {optimization.message}", result)
        return result


def solve_fluid_inclusion(problem: FluidInclusionProblem, **kwargs: object) -> FluidInclusionResult:
    """Convenience entry point for the stable public API."""
    return FluidInclusionSolver().solve(problem, **kwargs)
