import numpy as np
import random


def axisUnitVec(axis=1):
    if axis == 1:
        return np.array([1, 0, 0])
    if axis == 2:
        return np.array([0, 1, 0])
    if axis == 3:
        return np.array([0, 0, 1])
    else:
        print("Bruh, it's a 3-vec.")


def randomUnitVec():
    theta = random.random()
    r = random.random()
    x = np.sin(theta) * r
    y = np.cos(theta) * r
    z = -np.sqrt(abs(x**2 + y**2 - 1))
    return np.array([x, y, z])


def randomUnitVec4d(x0=1):
    """
    WARNING: This isn't technically a unit vector. The temporal component will blow the pythagorean relationship up, but that relationship along our temporal dimension doesn't appear to be trignomeric anyways...

    The temporal component will always have a magnitude of 1.
    """
    theta = random.random()
    r = random.random()
    x = np.sin(theta) * r
    y = np.cos(theta) * r
    z = -np.sqrt(abs(x**2 + y**2 - 1))
    return np.array([x0, x, y, z])
