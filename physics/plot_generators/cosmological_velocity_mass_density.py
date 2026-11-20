"""Plot derived cosmological velocity against mass density for all bodies."""

import matplotlib.pyplot as plt
import scipy.constants as spc

from physics.data.subsets import getDataFrame

offsetMap = {
    "Saturn": (-15, -10),
    "Uranus": (15, -20),
    "Jupiter": (15, -15),
    "Sun": (15, -10),
    "Neptune": (-40, 0),
    "Haumea": (-50, 0),
    "Pluto": (-35, 5),
    "Makemake": (0, -15),
    "Ceres": (10, -10),
    "Eris": (15, -10),
    "Moon": (-15, -15),
    "Mars": (-15, -15),
    "Venus": (-15, -15),
    "Mercury": (-35, -35),
    "Earth": (-30, -30),
}


def plotCosmologicalVelocityMassDensity():
    """Plot the derived cosmological velocity for every body in the dataset.

    The velocity is calculated as ``c * sqrt(1 - 1 / (1 + g / R)**2)``,
    where ``g`` is surface gravitational acceleration and ``R`` is radius.
    The dataframe is loaded without excluding outliers. Returns ``(fig, ax)``.
    """
    bodies = getDataFrame(drop_outliers=True)
    mass_density = bodies["mass_density"]
    velocity = (
        spc.c
        * (1 - 1 / (1 + bodies["g_acceleration"] / bodies["radius_m"]) ** 2) ** 0.5
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.scatter(
        mass_density,
        velocity,
        s=70,
        alpha=0.85,
        color="#0072B2",
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )

    for name, density, speed in zip(bodies["name"], mass_density, velocity):
        offset_x = offsetMap[name][0]
        offset_y = offsetMap[name][1]
        ax.annotate(
            name,
            (density, speed),
            xytext=(offset_x, offset_y),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(r"Mass density (kg m$^{-3}$)", fontsize=11)
    ax.set_ylabel(
        r"Derived cosmological velocity $c\sqrt{1 - 1/(1 + g/R)^2}$ (m s$^{-1}$)",
        fontsize=11,
    )
    ax.set_title("Derived Cosmological Velocity vs. Mass Density", fontsize=13)
    ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="major", color="#b0b0b0", linewidth=0.7, alpha=0.35)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig, ax


if __name__ == "__main__":
    fig, ax = plotCosmologicalVelocityMassDensity()
    fig.show()
    fig.waitforbuttonpress()
