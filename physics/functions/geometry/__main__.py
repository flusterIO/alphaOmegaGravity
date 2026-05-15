import numpy as np
import numpy.typing as npt

def hypotonuse(arr: npt.NDArray):
    return np.pow(np.pow(arr, 2).sum(), 0.5)

