"""Private formula-keyed facade used by the fluid-inclusion solver."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

from ._core import PitzerProperties
from ._database import ParameterRepository


@dataclass(frozen=True)
class SolutionProperties:
    temperature_k: float
    molalities: Mapping[str, float]
    ionic_strength: float
    osmotic_coefficient: float
    water_activity: float
    ln_activity_coefficients: Mapping[str, float]
    mean_activity_coefficients: Mapping[tuple[str, str], float]

    def activity_coefficient(self, formula: str) -> float:
        import math
        return math.exp(self.ln_activity_coefficients[formula])


class PropertyCalculator:
    """Calculate solution properties using records resolved from a repository."""

    def __init__(self, repository: ParameterRepository | None = None):
        self.repository = repository or ParameterRepository()

    def calculate(self, molalities: Mapping[str, float], temperature_k: float = 298.15) -> SolutionProperties:
        model_input, ids = self.repository.pitzer_input(molalities, temperature_k)
        result = PitzerProperties(model_input).calculate()
        id_to_formula = {species_id: formula for formula, species_id in ids.items()}
        ln_gamma = MappingProxyType({id_to_formula[key]: value for key, value in result.ln_gamma.items()})
        mean_gamma = MappingProxyType({(id_to_formula[cation], id_to_formula[anion]): value for (cation, anion), value in result.mean_activity_coeff.items()})
        return SolutionProperties(
            temperature_k=float(temperature_k),
            molalities=MappingProxyType({name: float(value) for name, value in molalities.items()}),
            ionic_strength=result.ionic_strength,
            osmotic_coefficient=result.osmotic_coefficient,
            water_activity=result.water_activity,
            ln_activity_coefficients=ln_gamma,
            mean_activity_coefficients=mean_gamma,
        )
