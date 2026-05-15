import numpy as np
import astropy.constants.codata2010 as cd



def lp(hbar = cd.hbar, G = cd.G, c = cd.c):
    np.sqrt(hbar * G / c**3)
