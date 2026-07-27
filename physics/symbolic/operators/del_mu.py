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


def del_mu(gmn, gst, V, mu=1, nu=1):
    r"""
    Returns the operator from Sean Caroll's book.
    $$
    \nabla_{\nu} V^{\mu} = \partial_\mu V^{\nu} + \Gamma^{\nu}_{\mu \lambda}V^{\lambda}
    $$
    """
    d_gmn_dxs = smp.diff(
        vecComponentFromIndex(V, nu), coordSysComponentFromIndex(coordSys, s)
    )
