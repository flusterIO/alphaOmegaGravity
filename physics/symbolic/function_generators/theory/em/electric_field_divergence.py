from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.variables import c


def electricFieldDivergence():
    """
    Since
    $$
    \\nabla \\cdot E = c \\left(\\nabla \\cdot B \\right)
    $$
    We find:

    $$
    \\nabla \\cdot E = 2 c \\frac{GM}{R^3}
    $$
    """
    return c * delta()
