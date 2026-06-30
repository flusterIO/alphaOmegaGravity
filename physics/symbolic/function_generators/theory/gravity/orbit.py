import sympy as smp
from physics.symbolic.variables import G, M, R


def orbital_velocity_theory(r=R, m=M, _G=G):
    return smp.sqrt(_G * m * (_G * m + 2 * r**3)) / r**2
