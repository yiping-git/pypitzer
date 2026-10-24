"""PyPitzer: fluid-inclusion quantification with the Pitzer activity model."""

from .fluid_inclusion import (
    ConvergenceError,
    DEFAULT_ELEMENT_SPECIES,
    FluidInclusionProblem,
    FluidInclusionResult,
    FluidInclusionSolver,
    SolverDiagnostics,
    solve_fluid_inclusion,
)

__all__ = [
    "ConvergenceError",
    "DEFAULT_ELEMENT_SPECIES",
    "FluidInclusionProblem",
    "FluidInclusionResult",
    "FluidInclusionSolver",
    "SolverDiagnostics",
    "solve_fluid_inclusion",
]

__version__ = "2.0.0"
