import sympy as smp
from physics.symbolic.variables import c, M, R, alpha, G


def massDensityTraditional():
    return M / (4 / 3 * smp.pi * R**3)


def massDensityEquivalence():
    """
    Returns the derived mass density and the traditional mass density in a tuple, s:

    `[derivedValue, traditionalValue]`
    """
    a = alpha
    phi = smp.GoldenRatio
    a = phi / (G * c * (a * (2 - a)))
    b = massDensityTraditional()
    return [a, b]
