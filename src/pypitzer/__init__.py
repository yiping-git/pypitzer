"""PyPitzer: focused Pitzer properties and fluid-inclusion quantification."""

from .core import (
    BinaryParam,
    Ion,
    LambdaParam,
    MissingParameterError,
    OutOfRangeError,
    PitzerError,
    PitzerInput,
    PitzerProperties,
    PitzerResult,
    PsiParam,
    ThetaParam,
    ZetaParam,
)
from .database import (
    DATABASE_SHA256,
    DATABASE_VERSION,
    DatabaseError,
    ParameterNotFoundError,
    ParameterRepository,
    Reaction,
    ReactionNotFoundError,
)
from .fluid_inclusion import (
    ConvergenceError,
    DEFAULT_ELEMENT_SPECIES,
    FluidInclusionProblem,
    FluidInclusionResult,
    FluidInclusionSolver,
    SolverDiagnostics,
    solve_fluid_inclusion,
)
from .properties import PropertyCalculator, SolutionProperties

__all__ = [
    "BinaryParam",
    "ConvergenceError",
    "DEFAULT_ELEMENT_SPECIES",
    "DATABASE_SHA256",
    "DATABASE_VERSION",
    "DatabaseError",
    "FluidInclusionProblem",
    "FluidInclusionResult",
    "FluidInclusionSolver",
    "Ion",
    "LambdaParam",
    "MissingParameterError",
    "OutOfRangeError",
    "ParameterNotFoundError",
    "ParameterRepository",
    "PitzerError",
    "PitzerInput",
    "PitzerProperties",
    "PitzerResult",
    "PropertyCalculator",
    "PsiParam",
    "Reaction",
    "ReactionNotFoundError",
    "SolutionProperties",
    "SolverDiagnostics",
    "ThetaParam",
    "ZetaParam",
    "solve_fluid_inclusion",
]

__version__ = "2.0.0"
