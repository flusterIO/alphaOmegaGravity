import numpy as np
from physics.symbolic.variables import Ex, Ey, Ez, Bx, By, Bz


def F_sub_mu_nu():
    return np.array(
        [[0, Ex, Ey, Ez], [-Ex, 0, Bz, -By], [-Ey, -Bz, 0, Bx], [-Ez, By, Bx, 0]]
    )


def F_super_mu_nu():
    return np.array(
        [[0, -Ex, -Ey, -Ez], [Ex, 0, Bz, -By], [Ey, -Bz, 0, Bx], [Ez, By, Bx, 0]]
    )


def LLF(L, F, mu=0, nu=0, sigma=0, rho=0):
    r"""
    Returns a tuple representing the left and right sides of this equation:

    $$
    \left( E^{z} \right)^{\prime} = \left( F^{0 z} \right)^{\prime} = L_{x}^{0}L^{z}_z F^{xz}
    $$
    """
    a = F[mu][nu]
    b = L[mu][sigma] * L[nu][rho] * F[sigma][rho]
    return [a, b]
