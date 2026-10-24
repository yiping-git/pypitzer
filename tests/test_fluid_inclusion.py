import pytest

import pypitzer
from pypitzer import (
    FluidInclusionProblem,
    solve_fluid_inclusion,
    supported_last_melting_phases,
)


def test_public_api_is_fluid_inclusion_only():
    assert set(pypitzer.__all__) == {
        "ConvergenceError",
        "DEFAULT_ELEMENT_SPECIES",
        "FluidInclusionProblem",
        "FluidInclusionResult",
        "FluidInclusionSolver",
        "LastMeltingPhase",
        "SolverDiagnostics",
        "solve_fluid_inclusion",
        "supported_last_melting_phases",
    }
    assert not hasattr(pypitzer, "PropertyCalculator")
    assert not hasattr(pypitzer, "PitzerProperties")
    assert not hasattr(pypitzer, "ParameterRepository")


def test_published_k_na_kcl_example():
    result = solve_fluid_inclusion(
        FluidInclusionProblem(
            ratios={"K+": 2.0},
            last_melting_phase="KCl",
            temperature_c=25.0,
        )
    )

    assert result.diagnostics.converged
    assert result.molalities["Na+"] == pytest.approx(1.8452598062034364, rel=1e-8)
    assert result.molalities["K+"] == pytest.approx(3.690519612406873, rel=1e-8)
    assert result.molalities["Cl-"] == pytest.approx(5.535779418610021, rel=1e-8)
    assert abs(result.diagnostics.charge_balance_residual) < 1e-9
    assert abs(result.diagnostics.equilibrium_residual) < 1e-9


def test_standard_ratio_must_be_one():
    with pytest.raises(ValueError, match="standard ratio"):
        FluidInclusionProblem(
            ratios={"Na+": 2.0, "K+": 1.0},
            last_melting_phase="KCl",
            temperature_c=25.0,
        )


def test_ice_last_melting_constraint():
    result = solve_fluid_inclusion(
        FluidInclusionProblem(
            ratios={},
            last_melting_phase="ice",
            temperature_c=-10.0,
        )
    )
    assert result.diagnostics.converged
    assert result.molalities["Na+"] == pytest.approx(result.molalities["Cl-"])
    assert result.last_melting_phase == "ice"
    assert result.temperature_c == -10.0


def test_element_ratios_are_converted_to_species_constraints():
    problem = FluidInclusionProblem.from_element_ratios(
        {"K": 2.0},
        last_melting_phase="sylvite",
        temperature_c=25.0,
    )
    assert problem.ratios == {"K+": 2.0}
    assert problem.internal_standard == "Na+"


def test_supported_phases_are_discoverable_for_pysills():
    phases = {phase.phase_id: phase for phase in supported_last_melting_phases()}
    assert len(phases) == 46
    assert phases["KCl"].formula == "KCl"
    assert phases["calcite"].formula == "CaCO3"
    assert phases["NaCl-2H2O"].formula == "NaCl·2H2O"
    assert phases["KCl-MgCl2-6H2O"].formula == "KCl·MgCl2·6H2O"
    assert phases["KCl"].temperature_min_c == pytest.approx(-90.0)
    assert phases["KCl"].temperature_max_c == pytest.approx(300.0)
    assert {"calcite", "aragonite", "vaterite"} <= phases.keys()


def test_unknown_phase_has_clear_error():
    with pytest.raises(Exception, match="Unknown last-melting phase"):
        solve_fluid_inclusion(
            FluidInclusionProblem(
                ratios={"K+": 2.0},
                last_melting_phase="not-a-phase",
                temperature_c=25.0,
            )
        )
