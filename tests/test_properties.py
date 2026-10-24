import math

import pytest

from pypitzer import ParameterNotFoundError, PropertyCalculator


def test_nacl_property_regression_at_25_celsius():
    result = PropertyCalculator().calculate({"Na+": 1.0, "Cl-": 1.0}, 298.15)

    assert result.ionic_strength == pytest.approx(1.0)
    assert result.osmotic_coefficient == pytest.approx(0.9374537448961822, rel=1e-11)
    assert result.water_activity == pytest.approx(0.9667870930693, rel=1e-11)
    assert result.ln_activity_coefficients["Na+"] == pytest.approx(-0.41707260461259255, rel=1e-11)
    assert result.activity_coefficient("Na+") == pytest.approx(math.exp(-0.41707260461259255))


def test_required_binary_parameter_is_not_silently_zeroed():
    with pytest.raises(ParameterNotFoundError, match="No binary parameter"):
        PropertyCalculator().calculate({"Unobtainium+": 1.0, "Cl-": 1.0}, 298.15)
