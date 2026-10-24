# PyPitzer 2

PyPitzer is a focused thermodynamic package for quantifying fluid-inclusion
compositions from microthermometry and LA-ICP-MS internal-standard ratios. It
provides a Pitzer property calculator, a bundled read-only parameter database,
and a constrained fluid-inclusion solver.

Version 2 is a clean, incompatible redesign of the published version 1 package.
The original implementation remains available from the `v1-legacy` tag and
the `legacy-v1` branch.

## Installation

```bash
python -m pip install pypitzer
```

## Fluid-inclusion calculation

Ratios are molar aqueous-species/internal-standard ratios. Na is the default
internal standard and chloride is solved together with the absolute Na
molality.

```python
from pypitzer import FluidInclusionProblem, solve_fluid_inclusion

problem = FluidInclusionProblem(
    ratios={"K+": 2.0},
    equilibrium="KCl(s) = K+(aq) + Cl-(aq)",
    temperature_c=25.0,
    initial_molalities=(1.0, 1.0),
)

result = solve_fluid_inclusion(problem)
print(result.molalities)
print(result.diagnostics)
```

This published regression gives approximately 1.84526 mol/kg Na,
3.69052 mol/kg K, and 5.53578 mol/kg Cl. Fluid-inclusion temperatures retain
the published implementation's 273.16 Celsius-to-Kelvin offset for exact
reproducibility; direct property calculations accept Kelvin explicitly.

PySILLS-style elemental ratios can be converted with the typed constructor:

```python
problem = FluidInclusionProblem.from_element_ratios(
    {"K": 2.0},
    equilibrium="KCl(s) = K+(aq) + Cl-(aq)",
    temperature_c=25.0,
)
```

The built-in mapping covers the database-supported cations. Use
`species_by_element` when a different oxidation state or aqueous species is
required.

The solver imposes both electroneutrality and the last-melting-phase
equilibrium. Failure raises `ConvergenceError`; passing `raise_on_failure=False`
returns the result with full diagnostics instead.

## Direct property calculation

```python
from pypitzer import PropertyCalculator

properties = PropertyCalculator().calculate(
    {"Na+": 1.0, "Cl-": 1.0},
    temperature_k=298.15,
)
print(properties.water_activity)
print(properties.activity_coefficient("Na+"))
```

The SQLite database is bundled as a versioned package resource and is always
opened using SQLite read-only and query-only modes. Database generation and
editing tools are intentionally not part of the public package.

## Citation

Liu, Y., Wagner, T., & Fußwinkel, T. (2024). An integrated approach for
quantifying fluid inclusion data combining microthermometry, LA-ICP-MS, and
thermodynamic modeling. *Chemical Geology, 644*, 121863.
