from physics.symbolic.function_generators.theory.velocity.vbar import vbarCosmological
from physics.symbolic.variables import M


def earthMomentum():
    return M * vbarCosmological()
