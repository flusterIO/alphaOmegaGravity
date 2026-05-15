from physics.types.numberlike import NumberLike
import numpy as np
import astropy.constants.codata2010 as cd


def gamma(v: NumberLike, c: NumberLike = cd.c):
    return 1 / np.sqrt(1 - v**2 / c**2)


def beta(v: NumberLike, c: NumberLike = cd.c):
    return v**2 / c**2
