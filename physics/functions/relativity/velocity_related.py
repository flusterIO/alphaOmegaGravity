from physics.functions.relativity.gravity_related import gravityGammaEquivalent
from physics.types.numberlike import NumberLike
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
import numpy as np
import scipy.constants as spc


def gammaFromV(vel, c=spc.c):
    """
    Einstein's equation for gamma
    """
    return 1 / (1 - vel**2 / c**2) ** (1 / 2)


def vFromGamma(gamma: NumberLike):
    """
    $$
    v = c \\sqrt{1 - \\frac{1}{\\gamma^2}}
    $$
    """
    return cd.c * (1 - 1 / gamma**2) ** (1 / 2)


def vbarCosmological(
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    c: NumberLike = cd.c.value,
    G: NumberLike = cd.G.value,
):
    """
    $$
    \\bar{v} = c \\sqrt{1 - \\frac{1}{\\left(1 + \\frac{g}{R}\\right)^2}}
    $$
    """
    return vFromGamma(gravityGammaEquivalent(M, R, c, G))


def vbarKinematic():
    return cd.c * np.sqrt(
        1 - 1 / (1 + 1 / 2 * cd.G.value * ac.M_earth.value / ac.R_earth.value**3) ** 2
    )
