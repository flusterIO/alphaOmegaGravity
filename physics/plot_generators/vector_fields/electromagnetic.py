"""Interactive visualizations of electric and magnetic vector fields."""

import numpy as np
import sympy as smp

from physics.symbolic.variables import x, y, z


def plot_electric_and_magnetic_fields(
    electric_field,
    magnetic_field=None,
    *,
    x_range=(-1.0, 1.0),
    y_range=(-1.0, 1.0),
    z_range=(-1.0, 1.0),
    points=7,
    show=True,
):
    """Plot electric and magnetic vector fields together in interactive 3D.

    Parameters
    ----------
    electric_field : callable, sympy.MatrixBase, or 3-component sequence
        A callable of ``(x, y, z)`` returning three symbolic components, or
        symbolic components ``(E_x, E_y, E_z)`` expressed in the shared
        coordinates from :mod:`physics.symbolic.variables`.
    magnetic_field : callable, sympy.MatrixBase, or 3-component sequence, optional
        A callable of ``(x, y, z)`` returning ``(B_x, B_y, B_z)``, or symbolic
        components in the same coordinates. If
        omitted, the plot assumes a static electric field and sets ``B = 0``.
        An electric field alone does not determine a magnetic field in
        general, so provide this argument for time-dependent or sourced fields.
    x_range, y_range, z_range : tuple of float, optional
        Coordinate limits for the sampled plot region.
    points : int, optional
        Number of sample locations along each coordinate axis (minimum 2).
    show : bool, optional
        Whether to open the interactive Plotly figure before returning it.

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive 3D figure with separately toggleable electric and
        magnetic vector-field traces.

    Notes
    -----
    Vector lengths are normalized independently at each sample point so the
    arrows show direction clearly. Their physical magnitudes are available in
    the hover data; the two fields are not directly comparable by arrow length.
    """
    try:
        from plotly import graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "Interactive vector-field plots require Plotly. Install plotly to use this function."
        ) from exc

    if not isinstance(points, int) or points < 2:
        raise ValueError("points must be an integer greater than or equal to 2")

    def as_components(field, name):
        if callable(field):
            field = field(x, y, z)
        if isinstance(field, smp.MatrixBase):
            components = list(field)
        else:
            try:
                components = list(field)
            except TypeError as exc:
                raise TypeError(f"{name} must be a 3-component symbolic vector") from exc
        if len(components) != 3:
            raise ValueError(f"{name} must contain exactly three components")
        return components

    def sample(field, name, mesh):
        components = as_components(field, name)
        try:
            evaluate = smp.lambdify((x, y, z), components, modules="numpy")
            values = evaluate(*mesh)
            return tuple(np.broadcast_to(np.asarray(value, dtype=float), mesh[0].shape) for value in values)
        except (TypeError, ValueError, NameError) as exc:
            raise ValueError(
                f"Could not evaluate {name}; use expressions in shared x, y, and z variables"
            ) from exc

    axes = [
        np.linspace(*x_range, points),
        np.linspace(*y_range, points),
        np.linspace(*z_range, points),
    ]
    mesh = np.meshgrid(*axes, indexing="ij")
    electric = sample(electric_field, "electric_field", mesh)
    if magnetic_field is None:
        magnetic_field = (smp.S.Zero, smp.S.Zero, smp.S.Zero)
    magnetic = sample(magnetic_field, "magnetic_field", mesh)

    figure = go.Figure()
    for label, vector, color in (
        ("Electric field E", electric, "#e4572e"),
        ("Magnetic field B", magnetic, "#247ba0"),
    ):
        components = np.stack(vector, axis=-1)
        magnitude = np.linalg.norm(components, axis=-1)
        nonzero = magnitude > 0
        normalized = np.zeros_like(components, dtype=float)
        normalized[nonzero] = components[nonzero] / magnitude[nonzero, None]
        figure.add_trace(
            go.Cone(
                x=mesh[0].ravel(),
                y=mesh[1].ravel(),
                z=mesh[2].ravel(),
                u=normalized[..., 0].ravel(),
                v=normalized[..., 1].ravel(),
                w=normalized[..., 2].ravel(),
                customdata=magnitude.ravel(),
                hovertemplate=(
                    f"{label}<br>x=%{{x:.3g}}, y=%{{y:.3g}}, z=%{{z:.3g}}"
                    "<br>magnitude=%{customdata:.5g}<extra></extra>"
                ),
                name=label,
                showscale=False,
                colorscale=[[0, color], [1, color]],
                sizemode="absolute",
                sizeref=0.35,
                anchor="tail",
            )
        )

    figure.update_layout(
        title="Electric and Magnetic Vector Fields",
        scene={
            "xaxis_title": "x",
            "yaxis_title": "y",
            "zaxis_title": "z",
            "aspectmode": "cube",
        },
        legend={"title": "Field"},
    )
    if show:
        figure.show()
    return figure
