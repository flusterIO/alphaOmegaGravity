import numpy as np
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
from physics.types.numberlike import NumberLike


def orbital_velocity(M=ac.M_earth, R=ac.R_earth, G=cd.G):
    return np.sqrt(G * M / R)


def escape_velocity(M=ac.M_earth, R=ac.R_earth, G=cd.G):
    return np.sqrt(2 * G * M / R)


def gravitationalPotentialEnergy(M=ac.M_earth, R=ac.R_earth, G=cd.G):
    """
    This was converted from a negative quantity to reflect the change in direction of the acceleration vector in this model.

    > Notice that there is no dependence on the mass of the body who's potential is supposed to change, making this nothing more than a useful mathematical tool that falls out of the symmetries provided by the conservation laws.
    """
    return G * M / R

def centripitalAcceleration(v: NumberLike, R = ac.R_earth):
    """
    v²/R
    """
    return v**2 / R

def periodOfCircularOrbit(M: NumberLike, R = ac.R_earth, G  = cd.G):
    return 2 * np.pi * np.sqrt(R**3 / (G * M))

def centripitalForce(v: NumberLike, M: NumberLike, R = ac.R_earth):
    return M * v**2 / R
