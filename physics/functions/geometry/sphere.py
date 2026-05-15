import scipy.constants as spc


def sphereSurfaceArea(r=1):
    return 4 * spc.pi * r**2


def sphereVolume(r=1):
    return 4 / 3 * spc.pi * r**3
