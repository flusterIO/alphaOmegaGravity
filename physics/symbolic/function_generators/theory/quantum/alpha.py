from astropy.constants.iau2015 import M_sun
from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.function_generators.theory.gravity.omega import omega
import sympy as smp


def alphaEquivalence():
    """
    1 ± √ω

    This yields α and 2-α, satisfying the equation original equation:

    α = 1 / (2 − α) ∛(d/dr δ)
    """
    f = smp.sqrt(1 - omega())
    return [1 + f, 1 - f]


def alphaSolar():
    """
    (1/2 δ)^1/3
    """
    _n = (smp.Rational(1, 2) * delta()) ** smp.Rational(1 / 3)
    return [1 + _n, 1 - _n]
