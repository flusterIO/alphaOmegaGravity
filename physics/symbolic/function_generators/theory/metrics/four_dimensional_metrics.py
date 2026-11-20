import sympy as smp
from physics.symbolic.variables import vbar


def fourDimensionalMetricTensor():
    return smp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
