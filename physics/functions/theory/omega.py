import numpy as np
import scipy.constants as spc
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
from physics.functions.theory.delta import delta
from physics.types.numberlike import NumberLike


def omega(M: NumberLike = ac.M_earth, R: NumberLike = ac.R_earth, G: NumberLike = cd.G):
    return np.power(-delta(M, R, G), 1 / 3)
