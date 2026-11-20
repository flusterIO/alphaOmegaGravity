"""Symbolic values for general relativity."""

import sympy as smp
from physics.symbolic.variables import G, R, M, c, g, height, Rearth


def gravitational_time_dilation(
    mass=M, radius=R, gravitational_constant=G, speed_of_light=c
):
    """Return the Schwarzschild proper-time rate for a static clock.

    The returned ratio is ``dτ/dt`` for a clock at radius ``radius`` outside
    a spherical, non-rotating mass, relative to a clock at infinity.
    """
    return smp.sqrt(
        1 - 2 * gravitational_constant * mass / (radius * speed_of_light**2)
    )


def poundRebkaDilationByHeight():
    return g * height / c**2


def frequencyRatioByHeight(r1=R, r2=Rearth):
    """
    Returns the frequency at r1 divided by the frequency at r2.
    """
    return smp.sqrt(
        gravitational_time_dilation(radius=r1) / gravitational_time_dilation(radius=r2)
    )
