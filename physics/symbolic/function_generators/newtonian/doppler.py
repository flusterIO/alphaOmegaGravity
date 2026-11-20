from physics.functions.relativity.velocity_related import vbarCosmological
from physics.symbolic.variables import h, c, G


def newtonianDopplerShift(v=vbarCosmological(), c=c):
    return 1 + v / c
