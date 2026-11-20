import matplotlib.pyplot as plt
import scipy.constants as spc
import astropy.constants.iau2015 as ac
import numpy as np
from pathlib import Path

from physics.functions import averageOfMarginErrorTuple
from physics.data.constants.milky_way import (
    solarCircleEnclosedMass,
    solarObitalVelocity,
    galacticRadiusObservational,
)


# Data taken from 'The Circular Velocity Curve of the Milky Way from 5 to 25 kpc' A. Eilers et al.
observationalGalacticRadii = np.array(
    [
        1.62615209e20,
        1.77117893e20,
        1.92237713e20,
        2.07666101e20,
        2.22785921e20,
        2.41299987e20,
        2.52716994e20,
        2.70922492e20,
        2.86042312e20,
        3.01162132e20,
        3.16590520e20,
        3.31710340e20,
        3.47138728e20,
        3.62567116e20,
        3.77995504e20,
        3.93115324e20,
        4.08235144e20,
        4.23972100e20,
        4.39400488e20,
        4.54828876e20,
        4.69640128e20,
        4.85685651e20,
        5.01114039e20,
        5.16542427e20,
        5.32279383e20,
        5.47707771e20,
        5.62827591e20,
        5.78255979e20,
        5.93067231e20,
        6.08187051e20,
        6.25466846e20,
        6.41203801e20,
        6.55397918e20,
        6.72677713e20,
        6.83169017e20,
        7.01374514e20,
        7.30071316e20,
        7.65865176e20,
    ]
)

# The observed orbital velocity at each radius from the same Eilers paper.
observationalGalacticVelocities = (
    np.array(
        [
            226.83,
            230.8,
            231.2,
            229.88,
            229.61,
            229.91,
            228.86,
            226.5,
            226.2,
            225.94,
            225.68,
            224.73,
            224.02,
            223.86,
            222.23,
            220.77,
            220.92,
            217.47,
            217.31,
            217.6,
            217.07,
            217.38,
            216.14,
            212.52,
            216.41,
            213.7,
            207.89,
            209.6,
            206.45,
            201.91,
            199.84,
            198.14,
            195.3,
            213.67,
            176.97,
            193.11,
            176.63,
            198.42,
        ]
    )
    * spc.kilo
)


def derivedOrbitalVelocity(G=spc.G, M=ac.M_sun.value, R=spc.au):
    return np.sqrt(G * M * (G * M + 2 * R**3)) / R**2


def plotOrbitalVelocity():
    """Plot the Milky Way rotation curve against observations and local values."""
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    radius_kpc = observationalGalacticRadii / (spc.kilo * spc.parsec)
    velocity_kms = observationalGalacticVelocities / spc.kilo
    ax.plot(
        radius_kpc,
        [
            derivedOrbitalVelocity(G=spc.G, M=solarCircleEnclosedMass, R=r) / spc.kilo
            for r in observationalGalacticRadii
        ],
        label="Relativistic model",
        color="#0072B2",
        linewidth=2.2,
    )
    ax.plot(
        radius_kpc,
        velocity_kms,
        label="Eilers et al. (2019)",
        color="#222222",
        linewidth=0,
        marker="o",
        markersize=4.2,
        markeredgewidth=0,
        alpha=0.8,
        zorder=3,
    )
    ax.plot(
        radius_kpc,
        [
            (spc.G * solarCircleEnclosedMass / r) ** (1 / 2) / spc.kilo
            for r in observationalGalacticRadii
        ],
        label="Newtonian point-mass model",
        color="#D55E00",
        linestyle="--",
        linewidth=2,
    )
    # Mark the Sun's galactocentric radius and circular velocity.
    _galacticRadiusObservational = averageOfMarginErrorTuple(
        galacticRadiusObservational
    )
    _observedGalacticVelocity = solarObitalVelocity
    solar_radius_kpc = _galacticRadiusObservational / (spc.kilo * spc.parsec)
    solar_velocity_kms = _observedGalacticVelocity / spc.kilo
    ax.errorbar(
        solar_radius_kpc,
        solar_velocity_kms,
        label="Solar neighborhood",
        color="#009E73",
        marker="*",
        markersize=11,
        capsize=3,
        elinewidth=1.2,
        linestyle="none",
        zorder=4,
    )
    ax.set_xlabel("Galactocentric radius (kpc)", fontsize=11)
    ax.set_ylabel(r"Circular velocity (km s$^{-1}$)", fontsize=11)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, which="major", color="#b0b0b0", linewidth=0.7, alpha=0.35)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax


if __name__ == "__main__":
    fig, ax = plotOrbitalVelocity()
    fig.show()
    fig.waitforbuttonpress()
