import sympy as smp


def alphaWyler():
    return (9 / (8 * smp.pi**4)) * (
        (smp.pi**5 / (2**4 * smp.factorial(5))) ** (smp.Rational(1, 4))
    )
