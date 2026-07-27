from sympy.matrices import Matrix
from physics.symbolic.function_generators.theory.velocity.vbar import vbarCosmological
from physics.symbolic.variables import S


def vbarCosmologicalXVector():
    return Matrix([
        vbarCosmological(), 
        0, 
        0
    ])
