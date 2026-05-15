from physics.symbolic.function_generators.theory.gravity.delta import delta
import sympy as smp


def omega():
    return delta() ** smp.Rational(1, 3)
