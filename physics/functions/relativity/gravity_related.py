from physics.types.numberlike import NumberLike
import numpy as np
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
from astropy import units as u

from physics.values.unit_corrected_values import unitCorrected_G


def gravityGammaEquivalent(
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    c: NumberLike = cd.c.value,
    G: NumberLike = unitCorrected_G.value,
):
    """
    $$
    v = c \\sqrt{1 - 1 / (1 + g/R)^2}
    $$
    """
    return 1 + G * M / R**3


def gammaEquivTD(
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    c: NumberLike = cd.c.value,
    G: NumberLike = unitCorrected_G.value,
):
    """
    ## $\\gamma$ equivialent [T]ime [D]ependent

    $$
    v = c \\sqrt{1 - 1 / (1 + 1/2 g/R)^2}
    $$
    """
    return 1 + 0.5 * 1 / (1 + 1 / 2 * G * M / R**3) ** 2
