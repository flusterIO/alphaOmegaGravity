from astropy.constants.iau2015 import M_sun, R_sun
import numpy as np
import astropy.units as u
import astropy.constants.codata2010 as cd
from physics.functions.relativity.velocity_related import vbarCosmological
from physics.functions.theory.delta import delta
from physics.functions.theory.omega import omega
from physics.values.unit_utilities import metersReduced
import scipy.constants as spc


def Gamma(alpha=spc.alpha):
    """
    1 / (2 - α)
    """
    return 1 / (2 - alpha)


def alphaIsolated():
    """
    Returns α as it would appear in the `alpha` equation after it was solved, instead of relying on inserting the observational value of α.

    1 ± √(1 −ω)
    """
    _n = (1 - omega()) ** (1 / 2)
    return [1 + _n, 1 - _n]


def alphaEquivalent():
    """
    This returns the right hand side of:

    α = 1/(2 - α) −δ^1/3

    using the observational value for alpha
    """
    _n = 1 / (2 - spc.alpha) * delta() ** (1 / 3)
    return _n


def alphaSolar():
    """
    This returns the right hand side of:

    1 ± α = (1/2 δ)^1/3
    """
    _n = (1 / 2 * delta(M_sun.value, R_sun.value)) ** (1 / 3)
    return [1 * metersReduced + _n, 1 * metersReduced - _n]


def spinorApproximation():
    return 1 / (2 * u.m / u.s - spc.alpha * u.m / u.s)


def vbarAlphaEquivalence():
    """
    Returns
    $$
    (v_0 / c)^2 \\approx (\\alpha (2 - \\alpha))^3
    $$

    This demonstrates that our velocity _through_ space is linearly proportional to the dilation _of_ space,
    again demonstrating that time and cosmic inflation are the same process as observed from different reference frames,
    once this divergent geometry is accounted for.
    """
    a = vbarCosmological() ** 2 / cd.c**2
    b = (spc.alpha * (2 - spc.alpha)) ** 3
    return [a, b]
