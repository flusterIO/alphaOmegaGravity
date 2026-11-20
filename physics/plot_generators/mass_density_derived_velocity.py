from physics.data.subsets import getDataFrame
import matplotlib.pyplot as plt
import scipy.constants as spc
from physics.data.constants.earth import cmb_velocity_3_tuple

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


def plotMassDensityDerivedVelocity():
    bodies = getDataFrame(drop_outliers=True)
    massDensity = bodies["mass_density"].to_numpy()  # pyright: ignore
    names = bodies["name"]
    derivedVelocities = (
        spc.c
        * (1 - 1 / (1 + 0.5 * spc.G * bodies["mass_kg"] / bodies["radius_m"] ** 3) ** 2)
        ** (1 / 2)
    ).to_numpy() / cmb_velocity_3_tuple[
        1
    ]  # pyright: ignore

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.scatter(
        massDensity,
        derivedVelocities,
        s=70,
        alpha=0.85,
        color="#0072B2",
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )

    for name, density, velocity in zip(names, massDensity, derivedVelocities):
        offset_x = offsetMap[name][0]
        offset_y = offsetMap[name][1]
        ax.annotate(
            name,
            (density, velocity),
            xytext=(offset_x, offset_y),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(r"Mass density (kg m$^{-3}$)", fontsize=11)
    ax.set_ylabel(
        r"Derived orbital velocity $\frac{v_{\text{derived}}}{v_{\text{observed}}}$",
        fontsize=11,
    )
    ax.set_title("Derived Orbital Velocity vs. Mass Density", fontsize=13)
    ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="major", color="#b0b0b0", linewidth=0.7, alpha=0.35)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig, ax


if __name__ == "__main__":
    fig, ax = plotMassDensityDerivedVelocity()
    fig.show()
    fig.waitforbuttonpress()
