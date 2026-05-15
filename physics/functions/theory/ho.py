import numpy as np
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd
import scipy.constants as spc
from physics.types.numberlike import NumberLike
from physics.values.unit_corrected_values import unitCorrected_G


def hbar0Earth():
    """
    If this model is correct and gravity is just a macroscopic charge, the flux of space-time through the manifold should act like the Compton wave length, defining a unit of length where space becomes time.

    While undoubtedly controversial, there is nothing coontradictory about the fact that their might be an entire Universe within our planet's core if the speed of causality scales inversely proportional to this spatial density. It is only our human derived limitations that places a limit on the complexity of the very small. Small, after all, is completely relative. It's up to the geologists to figure out the rest.

    It is for this reason that $g$ acts as the Compton wavelength of a massive body, as per unit time the radius elongates by a factor of g, converting time into space in the truest physical sense.
    """
    return (
        unitCorrected_G.value
        * ac.M_earth.value**2
        / ac.R_earth.value**2
        * cd.c.value
    )


def h0Earth():
    return unitCorrected_G * ac.M_earth**2 / ac.R_earth**2 * cd.c * 2 * spc.pi


def hOfLambda(λ: NumberLike):
    """
    2πcλ
    """
    return 2 * spc.pi * cd.c * λ


def tau_proton_mass_ratio(
    G=cd.G.value, M=ac.M_earth.value, R=ac.R_earth.value, c=cd.c.value
):
    return (2 * G * M / R**3 / (c * np.sqrt(1 - 1 / (1 + G * M / R**3) ** 2))) ** 3
