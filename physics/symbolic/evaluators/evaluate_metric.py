import sympy as smp

from physics.symbolic.variables import G, c


def check_einstein_equations(
    metric, coordinates, stress_energy=None, *, G=G, c=c, Lambda=0
):
    """Check G_mn + Lambda*g_mn = (8*pi*G/c**4)*T_mn.

    Returns the symbolic residual matrix and whether every component
    simplifies to zero. Use stress_energy=zero_matrix for vacuum.
    """
    g = smp.Matrix(metric)
    coords = tuple(coordinates)
    n = len(coords)

    if g.shape != (n, n):
        raise ValueError("Metric shape must match the number of coordinates.")
    if g != g.T:
        raise ValueError("Metric must be symmetric.")

    if stress_energy is None:
        stress_energy = smp.zeros(n)
    T = smp.Matrix(stress_energy)
    if T.shape != (n, n):
        raise ValueError("Stress-energy tensor shape must match the metric.")

    g_inv = g.inv()
    Gamma = [[[smp.S(0) for _ in range(n)] for _ in range(n)] for _ in range(n)]

    for a in range(n):
        for b in range(n):
            for d in range(n):
                Gamma[a][b][d] = smp.simplify(
                    sum(
                        g_inv[a, e]
                        * (
                            smp.diff(g[e, d], coords[b])
                            + smp.diff(g[e, b], coords[d])
                            - smp.diff(g[b, d], coords[e])
                        )
                        for e in range(n)
                    )
                    / 2
                )

    Ricci = smp.zeros(n)
    for b in range(n):
        for d in range(n):
            Ricci[b, d] = smp.simplify(
                sum(
                    smp.diff(Gamma[a][b][d], coords[a])
                    - smp.diff(Gamma[a][b][a], coords[d])
                    + sum(
                        Gamma[a][a][e] * Gamma[e][b][d]
                        - Gamma[a][d][e] * Gamma[e][b][a]
                        for e in range(n)
                    )
                    for a in range(n)
                )
            )

    R = smp.simplify(sum(g_inv[b, d] * Ricci[b, d] for b in range(n) for d in range(n)))
    Einstein = Ricci - g * R / 2

    if G is None or c is None:
        # Geometric-units form: G_mn + Lambda*g_mn = 8*pi*T_mn.
        coupling = 8 * smp.pi
    else:
        coupling = 8 * smp.pi * G / c**4

    residual = (Einstein + Lambda * g - coupling * T).applyfunc(smp.simplify)
    satisfied = all(smp.simplify(component) == 0 for component in residual)

    return {
        "satisfied": satisfied,
        "residual": residual,
        "ricci_tensor": Ricci,
        "ricci_scalar": R,
        "einstein_tensor": Einstein,
    }
