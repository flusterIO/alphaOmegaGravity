import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as spc
from physics.data.subsets import outliers



def plot_gravity_alpha_relation(
    csv_path="/Users/bigsexy/Desktop/notes/content/resources/csv/gravity_geometry_dataset.csv",
    show_labels=True,
    log_x=True,
    log_y=True,
    drop_rows=[],
):
    """
    Plot the proposed gravitational quantity

        omega = (2 G M / R^3)^(1/3)

    against compactness

        C = G M / (R c^2)

    and compare omega against several candidate alpha relationships.

    Parameters
    ----------
    csv_path : str
        Path to the gravity geometry CSV dataset.

    show_labels : bool
        Annotate individual bodies on the plot.

    log_x : bool
        Use logarithmic compactness axis.

    log_y : bool
        Use logarithmic omega axis.

    Returns
    -------
    pandas.DataFrame
        The dataframe containing the calculated quantities.

    matplotlib.figure.Figure
        The generated figure.
    """

    # ---------------------------------------------------------
    # Constants
    # ---------------------------------------------------------

    G = spc.G
    c = spc.c

    # Fine-structure constant
    alpha = spc.alpha

    # ---------------------------------------------------------
    # Load data
    # ---------------------------------------------------------

    df = pd.read_csv(csv_path)
    for r in drop_rows:
        df = df[(df["name"] != r)].copy()

    # Make sure the required columns exist.
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

    delta = 2.0 * G * M / R**3

    omega = np.cbrt(delta)

    compactness = G * M / (R * c**2)
    df["delta"] = delta
    df["omega"] = omega
    df["compactness"] = compactness
    relationships = {
        "alpha": alpha,
        "alpha / 2": alpha / 2.0,
        "alpha(2-alpha)": alpha * (2.0 - alpha),
        "1 / (2-alpha)": 1.0 / (2.0 - alpha),
    }

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.scatter(
        compactness,
        omega,
        s=70,
        alpha=0.85,
    )

    if show_labels:

        for i, row in df.iterrows():
            ax.annotate(
                row["name"],
                (
                    row["compactness"],
                    row["omega"],
                ),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=8,
            )

    colors = ["red", "green", "orange", "blue"]
    for i, (name, value) in enumerate(relationships.items()):
        ax.axhline(
            value,
            linestyle="--",
            linewidth=1.2,
            label=f"{name} = {value:.6g}",
            color=colors[i],
        )

    if log_x:
        ax.set_xscale("log")

    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(
        r"Compactness  $C = \frac{GM}{Rc^2}$",
        fontsize=12,
    )

    ax.set_ylabel(
        r"$\omega = \left(\frac{2GM}{R^3}\right)^{1/3}$",
        fontsize=12,
    )

    ax.set_title(
        "Gravitational Geometry vs. Compactness",
        fontsize=15,
    )
    ax.grid(
        True,
        which="both",
        alpha=0.25,
    )

    ax.legend(
        fontsize=9,
        loc="best",
    )

    fig.tight_layout()

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
    _, fig, ax = plot_gravity_alpha_relation(drop_rows=outliers)
    show_fig(fig)
