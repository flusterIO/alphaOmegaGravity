import numpy as np
import scipy.constants as spc
import astropy.constants.iau2015 as ac

from physics.functions.relativity.velocity_related import vbarCosmological
from physics.types.numberlike import OptionalFloat


def metricTensorTheory(
    x1: OptionalFloat = None, x2: OptionalFloat = None, x3: OptionalFloat = None
):
    """
    Note: This is intentionally only 3 dimensional. Time *is* the dilation of this vector in the aog model.
    Also, if you pass in one x, pass in all 3 otherwise the output won't make triginometric sense.
    
    This gives dilation per unit distance traveled by the observer.
    """
    defaultXDot = spc.G * ac.M_earth.value / ac.R_earth.value**2 / vbarCosmological()
    _x1 = x1 if x1 is not None else defaultXDot
    _x2 = x2 if x2 is not None else defaultXDot
    _x3 = x3 if x3 is not None else defaultXDot
    return np.array([[_x1, 0, 0], [0, _x2, 0], [0, 0, _x3]])

def metricTensorTheory4dAFrame(
    x0: OptionalFloat = None,
    x1: OptionalFloat = None,
    x2: OptionalFloat = None,
    x3: OptionalFloat = None,
):
    """
    This variant references the coordinate system of the observer in relative motion as observed from the observer at relative rest.

    Includes a time component with the theoretical metric tensor as set to be equivalent to gravity by default: 1 + g/R.
    Also, if you pass in one x, pass in all 3 otherwise the output won't make triginometric sense.

    This gives diilation per unit distance traveled by the observer.
    """
    defaultXDot = spc.G * ac.M_earth.value / ac.R_earth.value**2 / vbarCosmological()
    _x0 = x0 if x0 is not None else defaultXDot
    _x1 = x1 if x1 is not None else 1 / defaultXDot
    _x2 = x2 if x2 is not None else 1/ defaultXDot
    _x3 = x3 if x3 is not None else 1 / defaultXDot
    return np.array([[_x0, 0, 0, 0], [0, _x1, 0, 0], [0, 0,  _x2, 0], [0, 0, 0, _x3]])

def metricTensorTheory4dBFrame(
    x0: OptionalFloat = None,
    x1: OptionalFloat = None,
    x2: OptionalFloat = None,
    x3: OptionalFloat = None,
):
    """
    This variant references the coordinate system of the observer in relative motion as observed from the observer that observer's perspective.

    Includes a time component with the theoretical metric tensor as set to be equivalent to gravity by default: 1 + g/R.
    Also, if you pass in one x, pass in all 3 otherwise the output won't make triginometric sense.
    """
    defaultXDot = 1 + spc.G * ac.M_earth.value / ac.R_earth.value**3
    _x0 = x0 if x0 is not None else defaultXDot
    _x1 = x1 if x1 is not None else 1
    _x2 = x2 if x2 is not None else 1
    _x3 = x3 if x3 is not None else 1
    return np.array([[_x0, 0, 0, 0], [0, _x1, 0, 0], [0, 0, _x2, 0], [0, 0, 0, _x3]])



