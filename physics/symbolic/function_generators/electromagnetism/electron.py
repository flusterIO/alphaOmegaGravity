import astropy.constants.codata2010 as cd
import sympy as smp
from physics.symbolic.variables import e, me, c,  eo, alpha


def electron_radius():
    """
    The classical formula for the electron radius.

    Interestingly, the following symmetry can be found:

    $$
    \\theta R_\\oplus = \\frac{P}{2 \\pi}
    $$

    Where $\\theta = \\sin{\\frac{1}{{R_e}}}^1$
    """
    return e**2 / ((me * c**2 * 8 * smp.pi * eo)) * 2


def gFactor(_alpha = alpha):
    return 2 + (_alpha / smp.pi) + (_alpha / smp.pi)**2 + (_alpha / smp.pi)**3 + (_alpha / smp.pi)**4
