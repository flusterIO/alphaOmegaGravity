from physics.types.numberlike import NumberLike
import scipy.constants as spc
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
import numpy as np
import astropy.units as u


def orbital_velocity_theory(M=ac.M_earth.value, R=ac.R_earth.value, G=cd.G.value):
    return np.sqrt(G * M * (G * M + 2 * R**3)) / R**2 * u.m / u.s


def dt_from_orbital_vel(M=ac.M_earth.value, R=ac.R_earth.value, G=cd.G.value):
    return R**3 / (0.5 * G * M - R**3)


def cOfV(v: NumberLike):
    """
    Since
    $$
    \\alpha = \\frac{v_{\\electron}}{c}
    $$

    It's only logical to wonder what 'c' is for an orbiting body to be used in similar calculations,
    so this returns $c$ or $c \\alpha$
    """
    v / spc.alpha


