import sympy as smp
from physics.symbolic.variables import R


def sphereSurfaceArea(r=R):
    return 4 * smp.pi * r**2


def sphereVolume(r=R):
    return smp.Rational(4, 3) * smp.pi * r**3
