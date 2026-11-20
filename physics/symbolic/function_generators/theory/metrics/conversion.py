import sympy as smp


def metric_from_ds2(ds2, differentials):
    """
    Extract g_ij from a quadratic line element ds2.

    Parameters
    ----------
    ds2 : sympy expression
        Line element expressed as a quadratic polynomial
        in the differential symbols.
    differentials : sequence of sympy symbols
        Differential-coordinate symbols, in metric index order.

    Returns
    -------
    sympy.Matrix
        Metric tensor g_ij, satisfying ds2 = dx.T * g * dx.
    """
    differentials = tuple(differentials)

    g = smp.hessian(ds2, differentials) / 2
    g = g.applyfunc(smp.simplify)

    # Verify that ds2 is represented by the extracted quadratic form.
    d = smp.Matrix(differentials)
    reconstructed = smp.expand((d.T * g * d)[0])

    if smp.simplify(ds2 - reconstructed) != 0:
        raise ValueError(
            "ds2 must be a homogeneous quadratic expression "
            "in the supplied differential symbols."
        )

    return g
