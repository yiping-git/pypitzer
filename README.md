# PyPitzer 2

PyPitzer is a focused thermodynamic package for quantifying fluid-inclusion
compositions from microthermometry and LA-ICP-MS internal-standard ratios. It
provides a constrained fluid-inclusion solver backed by a private Pitzer
runtime and a bundled read-only parameter database.

Version 2 is a clean, incompatible redesign of the published version 1 package.
The original implementation remains available from the `v1-legacy` tag and
the `legacy-v1` branch.

## Installation

```bash
python -m pip install pypitzer
```

## Fluid-inclusion calculation

`ratios_to_internal_standard` contains molar aqueous-species/internal-standard
ratios. `Na+` is the default internal standard, but another charged species may
be selected. The standard itself is omitted because its ratio is implicitly
one. Chloride is the default balancing ion.

```python
from pypitzer import FluidInclusionProblem, solve_fluid_inclusion

problem = FluidInclusionProblem(
    ratios_to_internal_standard={"K+": 2.0},  # K/Na = 2
    last_melting_phase="KCl",
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
reproducibility.

PySILLS-style elemental ratios can be converted with the typed constructor:

```python
problem = FluidInclusionProblem.from_element_ratios(
    element_ratios_to_internal_standard={"K": 2.0},  # K/Na = 2
    last_melting_phase="KCl",
    temperature_c=25.0,
)
```

The built-in mapping covers the database-supported cations. Use
`species_by_element` when a different oxidation state or aqueous species is
required.

PySILLS can populate its phase selector directly from the discovery API:

```python
from pypitzer import supported_last_melting_phases

for phase in supported_last_melting_phases():
    print(phase.phase_id, phase.formula)
```

Only the stable `phase_id` is passed back to `FluidInclusionProblem`; chemical
reactions remain private implementation details. Use `formula` as the label in
the PySILLS interface; `phase_id` remains distinct where one formula has
multiple polymorphs.

### Supported last-melting phases

| Phase ID | Formula | Valid °C |
|---|---|---:|
| `2KCl-FeCl2-2H2O` | 2KCl·FeCl2·2H2O | -90 to 25 |
| `aragonite` | CaCO3 | -90 to 25 |
| `CaCl2-2MgCl2-12H2O` | CaCl2·2MgCl2·12H2O | -90 to 25 |
| `CaCl2-6H2O` | CaCl2·6H2O | -90 to 25 |
| `CaCO3-6H2O` | CaCO3·6H2O | -90 to 25 |
| `calcite` | CaCO3 | -90 to 25 |
| `CaMg(CO3)2` | CaMg(CO3)2 | -90 to 25 |
| `CaSO4` | CaSO4 | -90 to 25 |
| `FeCl2-4H2O` | FeCl2·4H2O | -90 to 25 |
| `FeCl2-6H2O` | FeCl2·6H2O | -90 to 25 |
| `FeCO3` | FeCO3 | -90 to 25 |
| `FeSO4-7H2O` | FeSO4·7H2O | -90 to 25 |
| `FeSO4-H2O` | FeSO4·H2O | -90 to 25 |
| `H2O` | H2O | -90 to 25 |
| `K2Mg(SO4)2-6H2O` | K2Mg(SO4)2·6H2O | -90 to 25 |
| `K2SO4` | K2SO4 | -90 to 25 |
| `KCl` | KCl | -90 to 300 |
| `KCl-MgCl2-6H2O` | KCl·MgCl2·6H2O | -90 to 25 |
| `KHCO3` | KHCO3 | -90 to 25 |
| `LiCl` | LiCl | -90 to 25 |
| `LiCl-2H2O` | LiCl·2H2O | -90 to 300 |
| `LiCl-3H2O` | LiCl·3H2O | -90 to 25 |
| `LiCl-5H2O` | LiCl·5H2O | -90 to 25 |
| `LiCl-H2O` | LiCl·H2O | -90 to 300 |
| `LiOH` | LiOH | 50 to 200 |
| `LiOH-H2O` | LiOH·H2O | 25 to 200 |
| `Mg(HCO3)(OH)-2H2O` | Mg(HCO3)(OH)·2H2O | -90 to 25 |
| `Mg5(CO3)4(OH)2-4H2O` | Mg5(CO3)4(OH)2·4H2O | -90 to 25 |
| `MgCl2-12H2O` | MgCl2·12H2O | -90 to 25 |
| `MgCl2-6H2O` | MgCl2·6H2O | -90 to 25 |
| `MgCl2-8H2O` | MgCl2·8H2O | -90 to 25 |
| `MgCO3` | MgCO3 | -90 to 25 |
| `MgCO3-5H2O` | MgCO3·5H2O | -90 to 25 |
| `MgSO4-12H2O` | MgSO4·12H2O | -90 to 25 |
| `MgSO4-6H2O` | MgSO4·6H2O | -90 to 25 |
| `MgSO4-7H2O` | MgSO4·7H2O | -90 to 25 |
| `Na2CO3-10H2O` | Na2CO3·10H2O | -90 to 25 |
| `Na2CO3-7H2O` | Na2CO3·7H2O | -90 to 25 |
| `Na2CO3-NaHCO3-2H2O` | Na2CO3·NaHCO3·2H2O | -90 to 25 |
| `Na2Mg(SO4)2-4H2O` | Na2Mg(SO4)2·4H2O | -90 to 25 |
| `Na2SO4` | Na2SO4 | -90 to 25 |
| `Na2SO4-10H2O` | Na2SO4·10H2O | -90 to 25 |
| `NaCl` | NaCl | -90 to 300 |
| `NaCl-2H2O` | NaCl·2H2O | -90 to 25 |
| `NaHCO3` | NaHCO3 | -90 to 25 |
| `vaterite` | CaCO3 | -90 to 25 |

The solver imposes both electroneutrality and the last-melting-phase
equilibrium. Failure raises `ConvergenceError`; passing `raise_on_failure=False`
returns the result with full diagnostics instead.

The property engine and SQLite repository are private implementation details;
PyPitzer does not expose a general-purpose Pitzer property API. The database is
always opened using SQLite read-only and query-only modes. Database generation
and editing tools are intentionally not part of the package.

## Citation

Liu, Y., Wagner, T., & Fußwinkel, T. (2024). An integrated approach for
quantifying fluid inclusion data combining microthermometry, LA-ICP-MS, and
thermodynamic modeling. *Chemical Geology, 644*, 121863.
