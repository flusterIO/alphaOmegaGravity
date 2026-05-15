from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.function_generators.theory.gravity.gamma_equiv import (
    gravityGammaEquivalent,
)
from physics.symbolic.function_generators.theory.gravity.omega import omega
from physics.symbolic.variables import alpha, M, R, c, gamma, G
import sympy as smp


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
