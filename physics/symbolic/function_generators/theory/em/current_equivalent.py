from physics.symbolic.function_generators.theory.gravity.delta import deltaExpanded
from physics.symbolic.function_generators.theory.velocity.vbar import vbarCosmological
from physics.symbolic.variables import R, G, M, alpha
import sympy as smp


def currentEquivalent():
    """
    Returns:
    δ⋅v₀
    """
    return deltaExpanded() * vbarCosmological()


def newtonsConstantFromAlpha():
    return R**3 / M * ((1 / (1 - (alpha * (2 - alpha)) ** 3)) ** (1 / 2) - 1)
