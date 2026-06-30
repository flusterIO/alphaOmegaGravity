from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.function_generators.theory.gravity.gamma_equiv import (
    gravityGammaEquivalent,
)
from physics.symbolic.function_generators.theory.gravity.omega import omega
from physics.symbolic.variables import alpha, M, R, c, gamma, G, v
import sympy as smp


def gammaTraditional(_v = v, _c = c):
    return 1 / smp.sqrt(1 - _v**2 / _c**2)



def vFromGamma(_gamma=gamma):
    """
    $$
    v = c \\sqrt{1 - \\frac{1}{\\gamma^2}}
    $$
    """
    return c * (1 - 1 / _gamma**2) ** smp.Rational(1, 2)


def vbarCosmological():
    """
    $$
    \\bar{v} = c \\sqrt{1 - \\frac{1}{\\left(1 + \\frac{g}{R}\\right)^2}}
    $$
    """
    return vFromGamma(gravityGammaEquivalent())


def vbarKinematic() -> smp.Pow:
    return c * smp.sqrt(1 - 1 / (1 + smp.Rational(1, 2) * G * M / R**3) ** 2)


def vbarAlphaEquivalence():
    """
    Returns
    $$
    (v_0 / c)^2 \\approx (\\alpha (2 - \\alpha))^3
    $$

    This demonstrates that our velocity _through_ space is linearly proportional to the dilation _of_ space,
    again demonstrating that time and cosmic inflation are the same process as observed from different reference frames,
    once this divergent geometry is accounted for.
    """
    a = vbarCosmological() ** 2 / c**2
    b = (alpha * (2 - alpha)) ** 3
    return [a, b]
