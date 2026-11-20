from physics.symbolic.variables import h, c, G, hbar
import sympy as smp


def plankMass():
    """
    $$
    \\sqrt{\\frac{h}{2 \\pi} \\frac{c}{G}}
    $$
    """
    return smp.sqrt(hbar * c / G)


def plankTime():
    return smp.sqrt(hbar * G / c**5)


def plankLength():
    return c * plankTime()


def plankEnergy():
    return hbar / plankTime()
