from physics.symbolic.function_generators.geometry.coord_system_utils import (
    coordSysComponentFromIndex,
    radialDistanceFromOrigin,
    radialUnitVec,
    vecComponentFromIndex,
)
from physics.symbolic.tensors_kinda.cristoffel import cristoffel
from physics.symbolic.variables import S, G, M
import numpy as np
import sympy as smp
from sympy.tensor.array import Array


def del_mu(gmn, gst, V, mu=1, nu=1, coordSys=S):
    r"""
    Returns the operator from Sean Caroll's book.
    $$
    \nabla_{\nu} V^{\mu} = \partial_\mu V^{\nu} + \Gamma^{\nu}_{\mu \lambda}V^{\lambda}
    $$
    """
    d_Vnu_dx_mu = smp.diff(
        vecComponentFromIndex(V, nu), coordSysComponentFromIndex(coordSys, mu)
    )

    Gamma_mu1_nu_V_1 = cristoffel(gmn, gst, mu, 1, nu) * vecComponentFromIndex(V, 1)
    Gamma_mu2_nu_V_2 = cristoffel(gmn, gst, mu, 2, nu) * vecComponentFromIndex(V, 2)
    Gamma_mu3_nu_V_3 = cristoffel(gmn, gst, mu, 3, nu) * vecComponentFromIndex(V, 3)

    return d_Vnu_dx_mu + Gamma_mu1_nu_V_1 + Gamma_mu2_nu_V_2 + Gamma_mu3_nu_V_3
