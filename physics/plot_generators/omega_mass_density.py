import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as spc
from physics.data.subsets import outliers


offsetMap = {
    "Saturn": 5,
    "Uranus": 5,
    "Jupiter": 12,
    "Sun": 18,
    "Neptune": 5,
    "Haumea": 10,
    "Pluto": 20,
    "Makemake": -14,
    "Ceres": 6,
    "Eris": 16,
    "Moon": 12,
    "Mars": -14,
    "Venus": -18,
    "Mercury": 10,
    "Earth": 20,
}


def plot(
    csv_path="/Users/bigsexy/Desktop/notes/content/resources/csv/gravity_geometry_dataset.csv",
    show_labels=True,
    log_x=True,
    log_y=True,
    drop_rows=[],
):

    G = spc.G
    c = spc.c

    alpha = spc.alpha

    df = pd.read_csv(csv_path)

    for r in drop_rows:
        df = df[(df["name"] != r)].copy()

    required = {
        "name",
        "mass_kg",
        "radius_m",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Dataset is missing required columns: {sorted(missing)}")

    # Remove objects without usable mass/radius.
    df = df.dropna(
        subset=[
            "mass_kg",
            "radius_m",
        ]
    ).copy()

    df = df[(df["mass_kg"] > 0) & (df["radius_m"] > 0)].copy()
    M = df["mass_kg"].to_numpy()
    R = df["radius_m"].to_numpy()

    mass_density = M / (4 / 3 * spc.pi * R**3)

    delta = 2.0 * G * M / R**3

    omega = np.cbrt(delta)

    compactness = G * M / (R * c**2)
    df["delta"] = delta
    df["omega"] = omega
    df["compactness"] = compactness
    df["mass_density"] = mass_density
    relationships = {
        r"$\alpha$": alpha,
        r"$\frac{\alpha}{2}$": alpha / 2.0,
        r"$\alpha(2-\alpha)$": alpha * (2.0 - alpha),
        r"$\frac{1}{(2-\alpha)}$": 1.0 / (2.0 - alpha),
    }

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)

    ax.scatter(
        mass_density,
        omega,
        s=70,
        alpha=0.85,
        color="#0072B2",
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )

    df = df.sort_values(by="mass_density", ignore_index=True)

    if show_labels:
        for i, row in df.iterrows():
            n = offsetMap[row["name"]] if row["name"] in offsetMap else 5
            ax.annotate(
                row["name"],
                (
                    row["mass_density"],
                    row["omega"],
                ),
                xytext=(-20, n),
                textcoords="offset points",
                fontsize=8,
            )

    colors = ["#D55E00", "#009E73", "#E69F00", "#0072B2"]
    for i, (name, value) in enumerate(relationships.items()):
        ax.axhline(
            value,
            linestyle="--",
            linewidth=1.5,
            label=f"{name} = {value:.6g}",
            color=colors[i],
        )

    if log_x:
        ax.set_xscale("log")

    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(
        r"Mass Density  $\rho = \frac{M}{\frac{4}{3}\pi R^3}$",
        fontsize=11,
    )

    ax.set_ylabel(
        r"$\omega = \left(\frac{2GM}{R^3}\right)^{1/3}$",
        fontsize=11,
    )

    ax.set_title(
        r"$\omega$ vs. Mass Density",
        fontsize=13,
    )
    ax.grid(True, which="major", color="#b0b0b0", linewidth=0.7, alpha=0.35)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        fontsize=9,
        loc="best",
        frameon=False,
    )

    # ---------------------------------------------------------
    # Print useful numerical results
    # ---------------------------------------------------------

    print()
    print("=" * 70)
    print("GRAVITATIONAL GEOMETRY SWEEP")
    print("=" * 70)
    print()

    print(
        df[
            [
                "name",
                "mass_kg",
                "radius_m",
                "compactness",
                "delta",
                "omega",
            ]
        ].to_string(index=False)
    )

    print()
    print("=" * 70)
    print("CANDIDATE ALPHA RELATIONSHIPS")
    print("=" * 70)
    print()

    for name, value in relationships.items():
        print(f"{name:20s} = {value:.12g}")

    print()
    print("=" * 70)

    return df, fig, ax


def show_fig(fig):
    fig.show()
    fig.waitforbuttonpress()


if __name__ == "__main__":
    _, fig, ax = plot(drop_rows=outliers)
    show_fig(fig)
