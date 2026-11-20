import sympy as smp
from physics.symbolic.variables import beta, v, c


def gammaExpanded(v=v):
    return 1 / smp.sqrt(1 - v**2 / c**2)


def relativisticDopplerShift():
    return smp.sqrt((1 - beta) / 1 + beta)
