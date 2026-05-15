import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
import scipy.constants as spc
import numpy as np


def plank_mass(h = spc.h, c = spc.c):
    return np.sqrt(h / (2 * spc.pi) * c / spc.G)
