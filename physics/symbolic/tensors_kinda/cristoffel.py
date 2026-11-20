import sympy as smp
from sympy.physics.vector import ReferenceFrame
from sympy.vector import CoordSys3D
from physics.symbolic.tensors_kinda import gmn
from physics.symbolic.variables import G, S
from physics.symbolic.function_generators.geometry.coord_system_utils import (
    coordSysComponentFromIndex,
)


def cristoffel(gmn, gst, m=1, n=1, t=1, s=1, coordSys: ReferenceFrame = S):
    r"""
    Returns
    $$
    \Gamma_{mn}^{t}
    $$

    Note: gst (upstairs) is the inverse of gmn (downstairs).
    """
    d_gsm_dxn = smp.diff(gmn[s - 1, t - 1], coordSysComponentFromIndex(coordSys, n))
    d_gsn_dxm = smp.diff(gmn[s - 1, n - 1], coordSysComponentFromIndex(coordSys, m))
    d_gmn_dxs = smp.diff(gmn[m - 1, n - 1], coordSysComponentFromIndex(coordSys, s))
    return smp.Rational(1, 2) * gst[s - 1, t - 1] * (d_gsm_dxn + d_gsn_dxm - d_gmn_dxs)


def cristoffelTensor(gmn, gst, coordSys: ReferenceFrame = S):
    r"""
    Returns
    $$
    \Gamma_{mn}^{t}
    $$

    Note: gst (upstairs) is the inverse of gmn (downstairs).
    """
    ms = []
    for m in range(3):
        ns = []
        for n in range(3):
            ts = []
            for t in range(3):
                ss = []
                for s in range(3):
                    d_gsm_dxn = smp.diff(
                        gmn[s, t], coordSysComponentFromIndex(coordSys, n + 1)
                    )
                    d_gsn_dxm = smp.diff(
                        gmn[s, n], coordSysComponentFromIndex(coordSys, m + 1)
                    )
                    d_gmn_dxs = smp.diff(
                        gmn[m, n], coordSysComponentFromIndex(coordSys, s + 1)
                    )
                    res = (
                        smp.Rational(1, 2)
                        * gst[s, t]
                        * (d_gsm_dxn + d_gsn_dxm - d_gmn_dxs)
                    )
                    ss.append(res)
                ts.append(ss)
            ns.append(ts)
        ms.append(ns)
    return ms
