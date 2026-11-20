"""Kaluza-Klein metric construction utilities."""

import sympy as smp
from physics.symbolic.variables import kappa


def _base_dimension(metric, name):
    """Return the supported base dimension, validating that metric is square."""
    if metric.rows != metric.cols or metric.rows not in (3, 4):
        raise ValueError(f"{name} must be a 3x3 or 4x4 matrix")
    return metric.rows


def _column_vector(vector, dimension, name):
    """Convert a vector to a column matrix and validate its dimension."""
    vector = smp.Matrix(vector)

    if vector.shape == (1, dimension):
        vector = vector.T
    if vector.shape != (dimension, 1):
        raise ValueError(f"{name} must contain {dimension} components")

    return vector


def kaluza_klein_tensor(metric_4d, potential_4d, coupling=kappa):
    r"""Return the Kaluza-Klein metric with one additional dimension.

    The convention is

    .. math::
       ds^2 = g_{\mu\nu} dx^\mu dx^\nu
              + (dy + \kappa A_\mu dx^\mu)^2,

    so the returned block matrix is
    ``[[g + kappa**2*A*A.T, kappa*A], [kappa*A.T, 1]]``.

    The base metric may be 3x3 or 4x4. The potential must have the
    corresponding number of components. The returned matrix is therefore
    4x4 or 5x5.
    """
    metric_4d = smp.Matrix(metric_4d)
    dimension = _base_dimension(metric_4d, "metric_4d")
    potential_4d = _column_vector(potential_4d, dimension, "potential_4d")

    coupled_potential = coupling * potential_4d
    metric_block = metric_4d + coupled_potential * coupled_potential.T
    return metric_block.row_join(coupled_potential).col_join(
        coupled_potential.T.row_join(smp.ones(1, 1))
    )


def derive_kk_metric(g, A, phi, kappa):
    """Construct the lifted Kaluza-Klein metric and its inverse.

    Parameters
    ----------
    g : 3x3 or 4x4 matrix
        Base metric.
    A : vector or matrix
        Covariant potential with one component per base dimension.
    phi : SymPy expression
        Scalar describing the additional dimension's scale.
    kappa : SymPy expression
        Electromagnetic coupling constant.

    Returns
    -------
    dict
        The lifted metric, its inverse, and its determinant. A 3D base
        produces 4x4 matrices; a 4D base produces 5x5 matrices.
    """
    g = smp.Matrix(g)
    dimension = _base_dimension(g, "g")
    A = _column_vector(A, dimension, "A")

    if g.equals(g.T) is not True:
        raise ValueError("The base metric must be symmetric.")

    # Kaluza-Klein metric components.
    Gbase = g + kappa**2 * phi**2 * A * A.T
    Gmix = kappa * phi**2 * A
    Gextra = smp.Matrix([[phi**2]])

    # Assemble the complete lifted metric.
    G = Gbase.row_join(Gmix).col_join(Gmix.T.row_join(Gextra))

    return {
        "metric": G.applyfunc(smp.simplify),
        "inverse_metric": G.inv().applyfunc(smp.simplify),
        "determinant": smp.factor(G.det()),
    }


def kk_metric_inputs_from_ds2(ds2, coords, extra, phi, kappa):
    """Extract arguments for ``derive_kk_metric`` from a line element.

    Parameters
    ----------
    ds2 : SymPy expression
        Expanded quadratic line element in the supplied differentials.
    coords : sequence
        Three or four base coordinates.
    extra : SymPy expression
        Coordinate for the additional dimension.
    phi : SymPy expression
        Scale field for the additional dimension.
    kappa : SymPy expression
        Electromagnetic coupling.

    Returns
    -------
    dict
        Arguments for ``derive_kk_metric(g, A, phi, kappa)``.
    """
    coords = tuple(coords)
    dimension = len(coords)
    if dimension not in (3, 4):
        raise ValueError("Exactly three or four base coordinates are required.")

    differentials = tuple(
        smp.Symbol(f"d_{coord}", commutative=True) for coord in (*coords, extra)
    )

    # Extract the lifted metric from the quadratic line element.
    G = (smp.hessian(ds2, differentials) / 2).applyfunc(smp.simplify)

    # Verify the line element is quadratic in the supplied differentials.
    differential_column = smp.Matrix(differentials)
    reconstructed = smp.expand((differential_column.T * G * differential_column)[0])
    if smp.simplify(ds2 - reconstructed) != 0:
        raise ValueError(
            "The line element is not a homogeneous quadratic form "
            "in the supplied differentials."
        )

    # Read the mixed and extra-dimensional blocks from the lifted metric.
    Gbase = G[:dimension, :dimension]
    Gmix = G[:dimension, dimension]
    Gextra = G[dimension, dimension]

    if smp.simplify(Gextra - phi**2) != 0:
        raise ValueError(
            "The extra-dimensional metric component does not equal phi**2."
        )

    # Invert the KK mixed-block relation to recover A, then remove its
    # contribution from the lifted base block to recover g.
    A = (Gmix / (kappa * phi**2)).applyfunc(smp.simplify)
    g = (Gbase - kappa**2 * phi**2 * A * A.T).applyfunc(smp.simplify)

    # Check that the recovered blocks reproduce the Kaluza-Klein ansatz.
    expected_Gbase = g + kappa**2 * phi**2 * A * A.T
    if not Gbase.equals(expected_Gbase):
        raise ValueError("The line element does not match the Kaluza-Klein ansatz.")

    if not g.equals(g.T):
        raise ValueError("The recovered base metric is not symmetric.")

    return {
        "g": g,
        "A": A,
        "phi": phi,
        "kappa": kappa,
    }
