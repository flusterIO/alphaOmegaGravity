from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.function_generators.theory.velocity.vbar import vbarCosmological


def currentEquivalent():
    """
    Returns:
    δ⋅v₀
    """
    return delta() * vbarCosmological()
