from sympy.vector import CoordSys3D
import sympy as smp

from physics.symbolic.variables import S


def vecComponentFromIndex(s, idx=1, coordSys=S):
    """
    Takes a **symbolic** vec.
    """
    if idx == 1:
        return s.components[coordSys.i]
    if idx == 2:
        return s.components[coordSys.j]
    if idx == 3:
        return s.components[coordSys.k]
    else:
        print(f"Invalid index! Found {idx}")


def coordSysComponentFromIndex(s: CoordSys3D, idx=1):
    if idx == 1:
        return s.x
    if idx == 2:
        return s.y
    if idx == 3:
        return s.z
    else:
        print(f"Invalid index! Found {idx}")


def radialDistanceFromOrigin(s: CoordSys3D):
    return (s.x**2 + s.y**2 + s.z**2) ** (smp.Rational(1, 2))


def radialUnitVec(s: CoordSys3D):
    mag = radialDistanceFromOrigin(s)
    return (s.x * s.i + s.y * s.j + s.z * s.k) / mag
