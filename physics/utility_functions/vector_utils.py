import numpy as np

def vectorMag(v):
    return np.sqrt(np.sum(v**2))


def vec4dsSquared(v):
    """
    Returns 
    ds^2 = -v_0^2 + v_1^2 + v_2^2 + v_3^2

    to match the negative temporal component of GR.
    """
    return -v[0]**2 + v[1]**2 + v[2]**2 + v[3]**2
