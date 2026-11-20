"""Symbolic force-ratio helpers involving photons."""

import sympy as smp

from physics.symbolic.variables import G, M, c, n_r, r


def photon_gravity_to_electromagnetic_force_ratio():
    r"""Return the ratio of radial photon accelerations from gravity and refraction.

    In the weak-field Newtonian approximation, use ``G*M/r**2`` for the
    inward gravitational acceleration. For an outward radial ray in an
    isotropic medium with refractive index ``n(r)``, its speed is ``c/n``;
    the magnitude of its radial coordinate acceleration is
    ``c**2*abs(dn/dr)/n**3``. The returned expression is their ratio.

    This compares coordinate accelerations in the stated approximations; it
    is not a vacuum electromagnetic force on a photon or a full GR result.
    """
    gravitational_acceleration = G * M / r**2
    electromagnetic_acceleration = (
        c**2 * smp.Abs(smp.diff(n_r, r)) / n_r**3
    )
    return gravitational_acceleration / electromagnetic_acceleration
