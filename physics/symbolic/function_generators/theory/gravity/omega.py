from physics.symbolic.function_generators.theory.gravity.delta import deltaExpanded
import sympy as smp


def omega():
    return deltaExpanded() ** smp.Rational(1, 3)
