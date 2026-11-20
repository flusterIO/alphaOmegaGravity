import scipy.constants as smpc
import sympy as smp

from physics.symbolic.variables import TDay, thetaInclination


def earth_rotation_matrix(dt, obliquity=thetaInclination):
    """
    Earth's rotation matrix in an ecliptic/inertial coordinate system.

    Parameters
    ----------
    dt : sympy expression or numeric
        Time interval in seconds.

    obliquity : sympy expression or numeric
        Earth's axial tilt relative to the ecliptic, in degrees.

    Returns
    -------
    sympy.Matrix
        3x3 rotation matrix for Earth's rotation over dt.
    """
    omega = 2 * smp.pi / TDay
    epsilon = obliquity

    theta = omega * dt

    # Earth's rotation axis in ecliptic coordinates
    axis = smp.Matrix(
        [
            0,
            smp.sin(epsilon),
            smp.cos(epsilon),
        ]
    )

    # Skew-symmetric cross-product matrix
    K = smp.Matrix(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ]
    )

    I = smp.eye(3)

    # Rodrigues' rotation formula
    return (
        I * smp.cos(theta) + (1 - smp.cos(theta)) * axis * axis.T + smp.sin(theta) * K
    )
