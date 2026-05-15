from physics.symbolic.variables import h, c, G
import sympy as smp

def plankMass():
    """
    $$
    \\sqrt{\\frac{h}{2 \\pi} \\frac{c}{G}}
    $$
    """
    return smp.sqrt(h / (2 * smp.pi) * c / G)

def plankLength():
    return smp.sqrt(h / (2 * smp.pi) * G / c**3)
