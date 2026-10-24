"""PyPitzer: fluid-inclusion quantification with the Pitzer activity model."""

from .fluid_inclusion import (
    ConvergenceError,
    DEFAULT_ELEMENT_SPECIES,
    FluidInclusionProblem,
    FluidInclusionResult,
    FluidInclusionSolver,
    LastMeltingPhase,
    SolverDiagnostics,
    solve_fluid_inclusion,
    supported_last_melting_phases,
)

__all__ = [
    "ConvergenceError",
    "DEFAULT_ELEMENT_SPECIES",
    "FluidInclusionProblem",
    "FluidInclusionResult",
    "FluidInclusionSolver",
    "LastMeltingPhase",
    "SolverDiagnostics",
    "solve_fluid_inclusion",
    "supported_last_melting_phases",
]

__version__ = "2.0.0"
