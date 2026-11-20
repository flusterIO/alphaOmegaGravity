import astropy.constants.codata2010 as cd
import sympy as smp
from physics.symbolic.variables import e, me, c,  eo, alpha, hbar, m_plank


def FeOverFg():
    return alpha * (m_plank / me)**2


def electron_radius():
    """
    The classical formula for the electron radius.

    Interestingly, the following symmetry can be found:

    $$
    \\theta R_\\oplus = \\frac{P}{2 \\pi}
    $$

    Where $\\theta = \\sin{\\frac{1}{{R_e}}}^1$
    """
    return alpha * hbar / (me * c)


def gFactor(_alpha = alpha):
    return 2 + (_alpha / smp.pi) + (_alpha / smp.pi)**2 + (_alpha / smp.pi)**3 + (_alpha / smp.pi)**4
