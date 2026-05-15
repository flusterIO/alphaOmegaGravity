import sympy as smp
from physics.symbolic.variables import e, hbar, mo, eo


def alphaTraditional():
    """
    The traditional representation of α
    """
    return 1 / (4 * smp.pi) * e**2 / hbar * (mo / eo) **smp.Rational(1, 2)
