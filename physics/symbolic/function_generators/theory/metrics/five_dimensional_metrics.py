import sympy as smp
from physics.symbolic.variables import R


def fiveDimensionalMetricTensor(metric_4d):
    r"""Extend a 4D metric by a compact circular fifth dimension.

    The product-metric ansatz is ``ds_5² = ds_4² + R² dχ²``, where ``χ`` is
    the dimensionless periodic coordinate on the compact circle and ``R`` is
    the shared symbolic compactification radius. The returned matrix is
    ``diag(metric_4d, R²)``. Gauge-field terms are not included because this
    function accepts only the 4D metric.
    """
    metric_4d = smp.Matrix(metric_4d)
    if metric_4d.shape != (4, 4):
        raise ValueError("metric_4d must be a 4x4 matrix")

    return metric_4d.row_join(smp.zeros(4, 1)).col_join(
        smp.zeros(1, 4).row_join(smp.Matrix([[R**2]]))
    )
