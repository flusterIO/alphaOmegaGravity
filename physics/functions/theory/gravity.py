import numpy as np
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
import scipy.constants as spc
import astropy.units as u
from physics.functions.statistics.statistic_utils import percentDiff
from physics.functions.theory.alpha import spinorApproximation
from physics.values import unit_corrected_values


def massDensityTraditional(M = ac.M_earth, R = ac.R_earth):
    return M / (4/3 * np.pi * R**3)

def massDensityEquivalence(G=unit_corrected_values.unitCorrected_G, c=cd.c):
    """
    Returns the derived mass density and the traditional mass density in a tuple, s:

    `[derivedValue, traditionalValue, percentError]`
    """
    _u = u.m / u.s
    a = spc.alpha * _u
    phi = spc.golden_ratio * u.m**3 / u.s**3
    a = (phi / (G * c * (a * (2 * _u - a))))
    b = massDensityTraditional()
    return [a, b, percentDiff(b.value, a.value)]

