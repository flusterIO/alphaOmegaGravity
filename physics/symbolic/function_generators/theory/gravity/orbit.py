import sympy as smp
from physics.symbolic.variables import G, M, R


def orbital_velocity_theory():
    return smp.sqrt(G * M * (G * M + 2 * R**3)) / R**2
