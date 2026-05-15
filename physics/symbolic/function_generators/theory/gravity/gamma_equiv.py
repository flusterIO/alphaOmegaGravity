from physics.symbolic.variables import M, R, G


def gravityGammaEquivalent():
    """
    $$
    v = c \\sqrt{1 - 1 / (1 + g/R)^2}
    $$
    """
    return 1 + G * M / R**3


def gammaEquivTD():
    """
    ## $\\gamma$ equivialent [T]ime [D]ependent

    $$
    v = c \\sqrt{1 - 1 / (1 + 1/2 g/R)^2}
    $$
    """
    return 1 + 0.5 * 1 / (1 + 1 / 2 * G * M / R**3) ** 2
