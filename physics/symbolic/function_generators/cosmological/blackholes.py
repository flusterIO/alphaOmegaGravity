from physics.symbolic.variables import G, eo, q, c, M
import sympy as smp


def Rs(m=M):
    return 2 * G * m / c**2


# def radiusQ():
#     return q**2 * G / (4 * smp.pi * eo * c**4)


# def chargedBlackHoleRadiusRatio():
#     """
#     Returns:

#     $$
#     \\frac{r_k}{r_s} = \\frac{1}{2} \\left(1 \\pm \\sqrt{1 - r_e^2} \\right)
#     $$

#     """
