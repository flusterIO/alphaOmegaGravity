from physics.symbolic.function_generators.quantum.alpha import alphaTraditional
from physics.symbolic.function_generators.theory.gravity.omega import omega
from physics.symbolic.variables import R
import sympy as smp


def magneticDivergence():
    """
    Returns the magnetic divergence equivalent:

    (α(2 − α))³
    """
    a = alphaTraditional()
    return (a * (2 - a)) ** 3


def magneticDivergenceIntegral():
    """
    Returns the radial integral of magnetic divergence:

    ∫ R³(α(2 − α))³ dR
    """
    a = alphaTraditional()
    return smp.integrate(R**3 * (a * (2 - a)) ** 3, R)
