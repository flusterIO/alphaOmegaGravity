import numpy as np
import astropy.constants.iau2015 as ac
from physics.symbolic.variables import G, M, R
import sympy as smp
import astropy.constants.codata2010 as cd
from physics.types.numberlike import NumberLike


def orbital_velocity():
    return np.sqrt(G * M / R)


def escape_velocity():
    return np.sqrt(2 * G * M / R)


def gravitationalPotentialEnergy():
    """
    This was converted from a negative quantity to reflect the change in direction of the acceleration vector in this model.

    > Notice that there is no dependence on the mass of the body who's potential is supposed to change, making this nothing more than a useful mathematical tool that falls out of the symmetries provided by the conservation laws.
    """
    return G * M / R

def centripitalAcceleration(v: smp.Symbol):
    """
    v²/R
    """
    return v**2 / R

def periodOfCircularOrbit():
    return 2 * np.pi * np.sqrt(R**3 / (G * M))

def centripitalForce(v: NumberLike):
    """
    Notice how there is no G in the centripital force equation.
    This alone should be indication enough that gravity is not a force, and that the net acceleration of the observer is zero.
    The math still works however as ds is equivalent to dx in magnitude and exactly opposite in direction, giving precisely the same change relative to the underlying coordinate system.
    """
    return M * v**2 / R
