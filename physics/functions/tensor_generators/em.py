from typing import Union
import numpy as np

from physics.functions.relativity.gamma_related import gamma
from physics.functions.relativity.velocity_related import vbarCosmological


def LorentzMatrix(v: Union[float, None] = None, _gamma: Union[float, None] = None):
    _v = v if v is not None else float(vbarCosmological().value)
    __gamma = _gamma if _gamma is not None else gamma(_v)
    return np.array(
        [
            [__gamma, -_v * __gamma, 0, 0],
            [-_v * __gamma, __gamma, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )


def ZTensor_T_A(v: Union[float, None] = None, _gamma: Union[float, None] = None):
    """
    Analogous to the Lorentz matrix, but giving the appropriate magnitude of T_B as
    measured in T_A.
    """
    _v = v if v is not None else float(vbarCosmological().value)
    __gamma = _gamma if _gamma is not None else gamma(_v)
    return np.array(
        [
            [__gamma, -_v / __gamma, 0, 0],
            [-_v / __gamma, __gamma, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )


def ZTensor_T_B(v: Union[float, None] = None, _gamma: Union[float, None] = None):
    """
    Analogous to the Lorentz matrix, but giving the appropriate magnitude of T_B as
    measured in T_B.
    """
    _v = v if v is not None else float(vbarCosmological().value)
    __gamma = _gamma if _gamma is not None else gamma(_v)
    return np.array(
        [
            [__gamma, -_v, 0, 0],
            [-_v, __gamma, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
