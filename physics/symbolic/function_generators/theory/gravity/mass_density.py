import sympy as smp
from physics.symbolic.variables import c, M, R, alpha, G, rho, dtheta, dt


def massDensityTraditional():
    return M / (smp.Rational(4, 3) * smp.pi * R**3)


def massDensityEquivalence():
    """
    Returns the derived mass density and the traditional mass density in a tuple, s:

    `[derivedValue, traditionalValue]`

    where the derived density is:
    Φ / (G c α (2 − α))
    """
    a = alpha
    phi = smp.GoldenRatio
    a = phi / (G * c * (a * (2 - a)))
    b = massDensityTraditional()
    return [a, b]


def GammaFromMassDensityEquivalence(mass_density=rho):
    return mass_density * alpha * G * c / smp.GoldenRatio


def GammaFromRotation():
    r"""
    Returns
    $$
    \sqrt{\frac{c \alpha^3}{\frac{d\theta}{dt} R}}
    $$
    which can be easily derived from the equation for Earth's rotation intersection
    """

    return smp.sqrt(c * alpha**3 / (dtheta / dt * R))
