import sympy as smp
from physics.symbolic.variables import G, M, R, c


def t_r_g():
    """
    Returns the velocity required to produce g when inserted into the vbar equation. This follows from the fact that Einstein's formula gives our escape velocity, which in turn is integrated radially to infinity.
    """
    return 1 / (1 - G**2 * M**2 / (R**2 * c**2))
