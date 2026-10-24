import pytest

from pypitzer import FluidInclusionProblem, solve_fluid_inclusion


def test_published_k_na_kcl_example():
    result = solve_fluid_inclusion(
        FluidInclusionProblem(
            ratios={"K+": 2.0},
            equilibrium="KCl(s) = K+(aq) + Cl-(aq)",
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
            equilibrium="KCl(s) = K+(aq) + Cl-(aq)",
            temperature_c=25.0,
        )


def test_ice_last_melting_constraint():
    result = solve_fluid_inclusion(
        FluidInclusionProblem(
            ratios={},
            equilibrium="H2O(s) = H2O(l)",
            temperature_c=-10.0,
        )
    )
    assert result.diagnostics.converged
    assert result.molalities["Na+"] == pytest.approx(result.molalities["Cl-"])
    assert result.properties.water_activity < 1.0


def test_element_ratios_are_converted_to_species_constraints():
    problem = FluidInclusionProblem.from_element_ratios(
        {"K": 2.0},
        equilibrium="KCl(s) = K+(aq) + Cl-(aq)",
        temperature_c=25.0,
    )
    assert problem.ratios == {"K+": 2.0}
    assert problem.internal_standard == "Na+"
