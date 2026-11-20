from sympy.physics.vector import ReferenceFrame
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


def coordSysComponentFromIndex(s: ReferenceFrame, idx=1):
    if idx == 1:
        return s[0]
    if idx == 2:
        return s[1]
    if idx == 3:
        return s[2]
    else:
        print(f"Invalid index! Found {idx}")


def radialDistanceFromOrigin(s: ReferenceFrame):
    return (s[0] ** 2 + s[1] ** 2 + s[2] ** 2) ** (smp.Rational(1, 2))


def radialUnitVec(s: CoordSys3D):
    mag = radialDistanceFromOrigin(s)
    return (s[0] * s.i + s[1] * s.j + s[2] * s.k) / mag
