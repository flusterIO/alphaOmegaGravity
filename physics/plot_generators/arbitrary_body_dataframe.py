"""Paper-styled plots for arbitrary pandas dataframes."""

from collections.abc import Mapping

import matplotlib.pyplot as plt
import pandas as pd


_PAPER_BLUE = "#0072B2"


def _plot_dataframe(
    dataframe,
    x,
    y,
    *,
    kind,
    label=None,
    xlabel=None,
    ylabel=None,
    title=None,
    color=_PAPER_BLUE,
    style=None,
    figsize=(7.2, 4.8),
    log_x=False,
    log_y=False,
    ax=None,
    **plot_properties,
):
    """Plot two dataframe columns with the project's paper figure styling.

    ``style`` may be a matplotlib format string (for example, ``"o"`` or
    ``"--"``) or a mapping of additional matplotlib artist properties.
    Remaining keyword arguments are passed to ``Axes.plot`` or
    ``Axes.scatter`` as appropriate.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError("dataframe must be a pandas DataFrame")
    if x not in dataframe.columns or y not in dataframe.columns:
        missing = [column for column in (x, y) if column not in dataframe.columns]
        raise ValueError(f"Dataframe is missing required columns: {missing}")
    if style is not None and not isinstance(style, (str, Mapping)):
        raise TypeError("style must be a matplotlib format string or mapping")

    x_values = dataframe[x]
    y_values = dataframe[y]
    valid = x_values.notna() & y_values.notna()
    x_values = x_values[valid]
    y_values = y_values[valid]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    properties = dict(plot_properties)
    if isinstance(style, Mapping):
        properties.update(style)

    properties.setdefault("color", color)
    if label is not None:
        properties.setdefault("label", label)

    if kind == "points":
        properties.setdefault("s", 70)
        properties.setdefault("alpha", 0.85)
        properties.setdefault("edgecolors", "white")
        properties.setdefault("linewidths", 0.6)
        properties.setdefault("zorder", 3)
        ax.scatter(x_values, y_values, **properties)
    else:
        properties.setdefault("linewidth", 2.2 if kind == "line" else 2.0)
        if kind == "dash":
            properties.setdefault("linestyle", "--")
        if isinstance(style, str):
            ax.plot(x_values, y_values, style, **properties)
        else:
            ax.plot(x_values, y_values, **properties)

    ax.set_xlabel(xlabel if xlabel is not None else str(x), fontsize=11)
    ax.set_ylabel(ylabel if ylabel is not None else str(y), fontsize=11)
    if title is not None:
        ax.set_title(title, fontsize=13)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    ax.grid(True, which="major", color="#b0b0b0", linewidth=0.7, alpha=0.35)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if label is not None:
        ax.legend(frameon=False, fontsize=9, loc="best")

    return fig, ax


def plot_dataframe_line(dataframe, x, y, **properties):
    """Plot one dataframe series as a solid line.

    Required inputs
    ---------------
    dataframe : pandas.DataFrame
        The table containing the values to plot. Each plotted series uses one
        x column and one y column; rows missing either value are skipped.
    x : str
        Name of the dataframe column to use for horizontal coordinates.
    y : str
        Name of the dataframe column to use for vertical coordinates.

    Common keyword properties
    -------------------------
    label : str, optional
        Name shown for this series in the legend. Supplying a label displays
        the legend. Give each overlaid series its own label.
    xlabel, ylabel : str, optional
        Text for the shared horizontal and vertical axis labels. Usually set
        these on the first series; later calls on the same axes can omit them.
        They default to the selected column names.
    title : str, optional
        Plot title.
    color : matplotlib color, optional
        Line color. Defaults to the paper blue ``"#0072B2"``.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw on. Use the axes returned by the first call to
        overlay more series; omit it to create a new figure.
    figsize : tuple, optional
        New figure size in inches, default ``(7.2, 4.8)``. Only used when
        ``ax`` is omitted.
    log_x, log_y : bool, optional
        Set the corresponding axis to logarithmic scale.
    style : str or mapping, optional
        A Matplotlib format string such as ``"-o"`` or a mapping of artist
        properties. Other Matplotlib ``Axes.plot`` options can be passed
        directly, such as ``linewidth=2``, ``linestyle=":"`` or
        ``marker="o"``.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes)
        The figure and axes, so later calls can use ``ax=ax``.

    Example
    -------
    >>> fig, ax = plot_dataframe_line(df, "time", "measurement",
    ...     label="Measurement", xlabel="Time (s)", ylabel="Position (m)")
    >>> plot_dataframe_line(df, "time", "model", label="Model",
    ...     color="#D55E00", ax=ax)
    """
    return _plot_dataframe(dataframe, x, y, kind="line", **properties)


def plot_dataframe_dash(dataframe, x, y, **properties):
    """Plot one dataframe series as a dashed line.

    ``dataframe``, ``x`` and ``y`` have the same meaning as in
    :func:`plot_dataframe_line`. Common keyword properties are ``label``
    (legend text), ``color``, ``xlabel``, ``ylabel``, ``title``, ``ax``,
    ``figsize``, ``log_x`` and ``log_y``. Pass ``ax`` from the first plot call
    to place this dash on the same figure. The default line style is ``"--"``;
    it can be changed with ``linestyle`` or ``style``. Matplotlib line options
    such as ``linewidth=1.5`` are also accepted.

    Returns ``(fig, ax)``. See :func:`plot_dataframe_line` for the full
    parameter descriptions and a multiple-series example.
    """
    return _plot_dataframe(dataframe, x, y, kind="dash", **properties)


def plot_dataframe_points(dataframe, x, y, **properties):
    """Plot dataframe values as points without connecting lines.

    ``dataframe``, ``x`` and ``y`` have the same meaning as in
    :func:`plot_dataframe_line`. Common keyword properties are ``label``
    (legend text), ``color``, ``xlabel``, ``ylabel``, ``title``, ``ax``,
    ``figsize``, ``log_x`` and ``log_y``. Use ``ax`` from a previous call to
    overlay points with lines or other point series. Matplotlib scatter options
    such as ``s=45`` (marker area), ``marker="s"``, and ``alpha=0.7`` can be
    passed directly.

    Returns ``(fig, ax)``. See :func:`plot_dataframe_line` for the full
    parameter descriptions.
    """
    return _plot_dataframe(dataframe, x, y, kind="points", **properties)


if __name__ == "__main__":
    # Replace these example values with a dataframe loaded from your data.
    example_data = pd.DataFrame(
        {
            "time": [0, 1, 2, 3, 4],
            "measured": [0.0, 0.9, 1.8, 2.7, 3.6],
            "model_a": [0.0, 1.0, 2.0, 3.0, 4.0],
            "model_b": [0.0, 0.7, 1.5, 2.4, 3.4],
        }
    )

    fig, ax = plot_dataframe_line(
        example_data,
        x="time",
        y="measured",
        label="Measurements",
        marker="o",
        linewidth=0,
    )
    plot_dataframe_line(
        example_data,
        x="time",
        y="model_a",
        label="Model A",
        color="#D55E00",
        ax=ax,
    )
    plot_dataframe_dash(
        example_data,
        x="time",
        y="model_b",
        label="Model B",
        color="#009E73",
        ax=ax,
    )

    fig.show()
