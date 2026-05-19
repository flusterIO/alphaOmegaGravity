from physics.symbolic.variables import alpha, au
import sympy as smp


def earthCEquivalentApproximation():
    """
    Since

    $$
    v = c \\alpha
    $$

    for an electron, and this model concludes that gravity is just
    the sum of charge as applied to space itself, this function gives
    the equivalent of $c$ in the function above for our Earth based parameters.
    """
    return 2 * smp.pi * au / (365.25 * 86400) / alpha
