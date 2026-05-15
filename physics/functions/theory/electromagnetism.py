import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
from physics.functions.relativity.velocity_related import vbarCosmological, vbarKinematic
from physics.functions.theory.alpha import spinorApproximation
from physics.functions.theory.delta import delta
from physics.types.numberlike import NumberLike
from physics.values import unit_corrected_values


def electricFieldDivergence(
    c = cd.c,
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    """
    Since
    $$
    \\nabla \\cdot E = c \\left(\\nabla \\cdot B \\right)
    $$
    We find:

    $$
    \\nabla \\cdot E = 2 c \\frac{GM}{R^3}
    $$
    """
    return c * delta(M, R, G)


def rotationalVelocityFromDivergence(
    c = cd.c,
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    return electricFieldDivergence(c, M, R, G) * spinorApproximation()


def earthChargeDensity(
    c = cd.c,
    mu0 = cd.mu0,
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    return electricFieldDivergence(c, M , R, G) * mu0

def ExB(
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    """
    Returns the magnitude of E⃗ × B⃗ where |B| = δR and |E| = δcR
    """
    _delta = delta(M, R, G)
    return _delta * ac.R_earth * _delta * cd.c * ac.R_earth

def poyntingMomentum(
    c = cd.c,
    v = vbarCosmological(),
    M: NumberLike = ac.M_earth.value,
    R: NumberLike = ac.R_earth.value,
    G: NumberLike = unit_corrected_values.unitCorrected_G.value,
):
    return ExB(M, R, G) / c**2 * v
