import sympy as smp
from physics.symbolic.variables import e, h, mo, eo


def alphaTraditional():
    """
    The traditional representation of α
    """
    return smp.Rational(1, 2) * e**2 / h * smp.sqrt(mo / eo)
