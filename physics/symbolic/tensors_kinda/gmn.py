from physics.symbolic.function_generators.geometry.coord_system_utils import (
    radialDistanceFromOrigin,
)
from physics.symbolic.function_generators.theory.velocity.vbar import vbarCosmological
from physics.symbolic.variables import S, G, M, R, Rearth
from sympy.tensor.array import Array


def g_mn():
    r = radialDistanceFromOrigin(S)
    dl = G * M / r**2 / vbarCosmological().subs([[R, Rearth]])
    return Array([[dl, 0, 0], [0, dl, 0], [0, 0, dl]])


def gSuper_mn():
    r = radialDistanceFromOrigin(S)
    one_over_dl = vbarCosmological().subs([[R, Rearth]]) / (G * M / r**2)
    return Array([[one_over_dl, 0, 0], [0, one_over_dl, 0], [0, 0, one_over_dl]])
