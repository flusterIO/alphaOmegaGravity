from physics.symbolic.variables import M, R, G, alpha, e, h, mo, eo
import sympy as smp


def _Gamma(_alpha=alpha):
    """
    1 / (2 - α)
    """
    return 1 / (2 - _alpha)


def gammaBar(r=R, m=M, _G=G):
    """
    $$
    v = c \\sqrt{1 - 1 / (1 + g/R)^2}
    $$
    """
    return 1 + _G * m / r**3


def gammaEquivTD():
    """
    ## $\\gamma$ equivialent [T]ime [D]ependent

    $$
    v = c \\sqrt{1 - 1 / (1 + 1/2 g/R)^2}
    $$
    """
    return 1 + 0.5 * 1 / (1 + 1 / 2 * G * M / R**3) ** 2


def gravitationalAccelFromElectromagnticParams(_alpha=alpha):
    return smp.Rational(1, 2) * (_alpha * (2 - _alpha)) ** 3
