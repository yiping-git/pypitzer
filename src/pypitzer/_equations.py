"""Numerical building blocks for the Pitzer model.

This module is private.  It contains only the equations needed by the runtime
property calculator and database coefficient evaluation.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np


def evaluate_coefficients(values: Sequence[float], temperature_k: float, equation: int) -> float:
    """Evaluate one of the coefficient forms used by the bundled database."""
    a = list(values)
    t = temperature_k
    if equation == 0:
        return a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] / t + a[5] * np.log(t)
    if equation == 1:
        return a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] / t + a[5] * np.log(t) + a[6] / t**2
    if equation == 2:
        return a[1] + a[2] * t + a[3] * t**2 + a[4] * t**3 + a[5] / t + a[6] * np.log(t) + a[7] / (t - 263.0) + a[8] / (680.0 - t)
    if equation == 3:
        return a[0] + a[1] * t + a[2] / t + a[3] * np.log(t) + a[4] / (t - 263.0) + a[5] * t**2 + a[6] / (680.0 - t) + a[7] / (t - 227.0)
    if equation == 4:
        tr = 298.15
        return a[1] + a[2] * (1.0 / t - 1.0 / tr) + a[3] * np.log(t / tr) + a[4] * (t - tr) + a[5] * (t**2 - tr**2) + a[6] * np.log(t - 260.0)
    if equation == 5:
        log10_k = a[0] + a[1] * t + a[2] / t + a[3] * np.log10(t) + a[4] / t**2
        return float(log10_k * np.log(10.0))
    if equation == 6:
        return a[0] + a[1] * t + a[2] / t + a[3] * np.log(t) + a[4] / t**2 + a[5] * t**2
    raise ValueError(f"Unsupported parameter equation {equation}")


def a_phi_spencer(temperature_k: float) -> float:
    values = [86.6836498, 0.0848795942, -8.88785150e-5, 4.88096393e-8, -1327.31477, -17.6460172]
    return float(evaluate_coefficients(values, temperature_k, 0))


def g_func(x: float) -> float:
    if x == 0.0:
        return 1.0
    return float(2.0 * (1.0 - (1.0 + x) * np.exp(-x)) / x**2)


def g_func_prime(x: float) -> float:
    if x == 0.0:
        return 0.0
    return float(-2.0 * (1.0 - (1.0 + x + x**2 / 2.0) * np.exp(-x)) / x**2)


def get_f_gamma(a_phi: float, ionic_strength: float) -> float:
    b = 1.2
    root_i = np.sqrt(ionic_strength)
    return float(-a_phi * (root_i / (1.0 + b * root_i) + (2.0 / b) * np.log(1.0 + b * root_i)))


_J_REGION_1 = np.array([1.925154015, -0.060076478, -0.029779077, -0.0072995, 0.000388261, 0.000636875, 3.65836e-05, -4.5037e-05, -4.5379e-06, 2.93771e-06, 3.96566e-07, -2.021e-07, -2.52678e-08, 1.35226e-08, 1.22941e-09, -8.21969e-10, -5.0847e-11, 4.6333e-11, 1.943e-12, -2.563e-12, -1.0991e-11])
_J_REGION_2 = np.array([0.628023321, 0.462762985, 0.150044637, -0.028796058, -0.036552746, -0.001668088, 0.00651984, 0.001130378, -0.000887171, -0.000242108, 8.72945e-05, 3.46821e-05, -4.58377e-06, -3.54868e-06, -2.50454e-07, 2.16992e-07, 8.07796e-08, 4.55856e-09, -6.94476e-09, -2.84926e-09, 2.37816e-10])


def _j_function(x: float) -> tuple[float, float]:
    if x <= 0.0:
        return 0.0, 0.25
    if x <= 1.0:
        z = 4.0 * x ** 0.2 - 2.0
        dz_dx = 0.8 * x ** -0.8
        coefficients = _J_REGION_1
    else:
        z = (40.0 / 9.0) * x ** -0.1 - 22.0 / 9.0
        dz_dx = -(4.0 / 9.0) * x ** -1.1
        coefficients = _J_REGION_2
    b = np.zeros(23)
    d = np.zeros(23)
    for k in range(20, -1, -1):
        b[k] = z * b[k + 1] - b[k + 2] + coefficients[k]
        d[k] = b[k + 1] + z * d[k + 1] - d[k + 2]
    return float(0.25 * x - 1.0 + 0.5 * (b[0] - b[2])), float(0.25 + 0.5 * dz_dx * (d[0] - d[2]))


def get_e_theta(z_m: int, z_n: int, a_phi: float, ionic_strength: float) -> tuple[float, float]:
    if z_m == z_n or ionic_strength == 0.0:
        return 0.0, 0.0
    def x(z1: int, z2: int) -> float:
        return 6.0 * z1 * z2 * a_phi * math.sqrt(ionic_strength)
    x_mn, x_mm, x_nn = x(z_m, z_n), x(z_m, z_m), x(z_n, z_n)
    j_mn, jp_mn = _j_function(x_mn)
    j_mm, jp_mm = _j_function(x_mm)
    j_nn, jp_nn = _j_function(x_nn)
    e_theta = (z_m * z_n / (4.0 * ionic_strength)) * (j_mn - 0.5 * j_mm - 0.5 * j_nn)
    derivative = -e_theta / ionic_strength + (z_m * z_n / (8.0 * ionic_strength**2)) * (x_mn * jp_mn - 0.5 * x_mm * jp_mm - 0.5 * x_nn * jp_nn)
    return float(e_theta), float(derivative)
