from physics.symbolic.function_generators.geometry.coord_system_utils import (
    coordSysComponentFromIndex,
    radialDistanceFromOrigin,
    radialUnitVec,
)
from physics.symbolic.tensors_kinda.cristoffel import cristoffel
from physics.symbolic.variables import S, G, M
import numpy as np
import sympy as smp
from sympy.tensor.array import Array


def curvatureTensor(gmn, gst, s=1, r=1, n=1, t=1, p=1, coordSys=S):
    r"""
    Returns
    $$
    R_{srn}^{t}
    $$
    """
    Gamma_sn_t = cristoffel(gmn, gst, s, n, t)
    d_Gamma_sn_t_dx_r = smp.diff(Gamma_sn_t, coordSysComponentFromIndex(coordSys, r))

    Gamma_rn_t = cristoffel(gmn, gst, r, n, t)
    d_Gamma_rn_t_dx_s = smp.diff(Gamma_rn_t, coordSysComponentFromIndex(coordSys, s))

    Gamma_sn_p = cristoffel(gmn, gst, s, n, p)
    Gamma_pr_t = cristoffel(gmn, gst, p, r, t)

    Gamma_rn_p = cristoffel(gmn, gst, r, n, p)
    Gamma_ps_t = cristoffel(gmn, gst, p, s, t)

    return d_Gamma_sn_t_dx_r - d_Gamma_rn_t_dx_s + Gamma_sn_p * Gamma_pr_t - Gamma_rn_p * Gamma_ps_t
