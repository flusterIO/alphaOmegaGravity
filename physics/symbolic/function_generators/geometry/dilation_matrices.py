import sympy as smp
from physics.symbolic.variables import phi, R, G, M, Rhat, dt, c, x, y, z


def dilationMatrixRadialProjection():
    return G * M / R**2 * Rhat


def dilation_matrix(r=None):
    if r is None:
        r = dilationMatrixRadialProjection()
    return smp.Matrix([[r, 0, 0], [0, r, 0], [0, 0, r]])


def radial_dilation_metric_tensor(dilation=None):
    r"""Return a Cartesian 4D metric with radial tangential-plane dilation.

    The default length scale is
    ``1/sqrt(1 - 2*G*M/(r*c**2))``, where ``r`` is computed from Cartesian
    coordinates ``(x, y, z)``. Radial separations are unchanged, while
    separations tangent to an equi-gravitational surface are scaled by this
    factor. The spatial metric is ``n*n.T + dilation**2*(I - n*n.T)``, where
    ``n`` is the Cartesian radial unit vector. The time component uses unit
    lapse and the ``(-, +, +, +)`` signature.
    """
    radius = smp.sqrt(x**2 + y**2 + z**2)
    if dilation is None:
        dilation = 1 / smp.sqrt(1 - 2 * G * M / (radius * c**2))

    radial_vector = smp.Matrix([x, y, z]) / radius
    radial_projector = radial_vector * radial_vector.T
    spatial_metric = radial_projector + dilation**2 * (
        smp.eye(3) - radial_projector
    )
    return smp.BlockMatrix(
        [[smp.Matrix([[-1]]), smp.zeros(1, 3)], [smp.zeros(3, 1), spatial_metric]]
    ).as_explicit()
