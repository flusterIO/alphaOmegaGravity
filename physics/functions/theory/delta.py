import astropy.constants.iau2015 as ac
import astropy.units as u
from physics.values.unit_utilities import metersPerSecond
from physics.types.numberlike import NumberLike
from physics.values import unit_corrected_values


def delta(
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    return 2 * G * M / R**3 * u.m
