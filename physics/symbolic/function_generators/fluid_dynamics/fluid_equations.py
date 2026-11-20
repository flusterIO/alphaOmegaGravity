"""Symbolic fluid-mechanics quantities and equations.

Vector components use the Cartesian order ``(x, y, z)``. The Navier–Stokes
helpers describe a Newtonian fluid with constant dynamic viscosity; the
momentum equation is written in its incompressible, constant-density form.
"""

import sympy as smp

from physics.symbolic.variables import (
    R,
    f_x,
    f_y,
    f_z,
    g,
    mu,
    p,
    rho,
    rho_field,
    t,
    u_x,
    u_y,
    u_z,
    ux,
    uy,
    x,
    y,
    z,
)

_COORDINATES = (x, y, z)
_DEFAULT_VELOCITY = (u_x, u_y, u_z)
_DEFAULT_BODY_FORCE = (f_x, f_y, f_z)


def velocity_vector(velocity=_DEFAULT_VELOCITY):
    """Return the Cartesian velocity as a SymPy column matrix."""
    return smp.Matrix(velocity)


def gradient(scalar, coordinates=_COORDINATES):
    """Return the spatial gradient of a scalar expression as a column vector."""
    return smp.Matrix([smp.diff(scalar, coordinate) for coordinate in coordinates])


def divergence(vector, coordinates=_COORDINATES):
    """Return the divergence of a three-component vector field."""
    components = tuple(vector)
    if len(components) != len(coordinates):
        raise ValueError("vector and coordinate dimensions must match")
    return sum(
        smp.diff(component, coordinate)
        for component, coordinate in zip(components, coordinates)
    )


def laplacian(scalar, coordinates=_COORDINATES):
    """Return the Cartesian scalar Laplacian, the sum of second spatial derivatives."""
    return sum(smp.diff(scalar, coordinate, 2) for coordinate in coordinates)


def vector_laplacian(vector, coordinates=_COORDINATES):
    """Apply the scalar Cartesian Laplacian independently to each vector component."""
    return smp.Matrix([laplacian(component, coordinates) for component in vector])


def material_derivative(
    scalar, velocity=_DEFAULT_VELOCITY, time=t, coordinates=_COORDINATES
):
    r"""Return the material derivative ``∂scalar/∂t + u·∇scalar``."""
    components = tuple(velocity)
    if len(components) != len(coordinates):
        raise ValueError("velocity and coordinate dimensions must match")
    return smp.diff(scalar, time) + sum(
        component * smp.diff(scalar, coordinate)
        for component, coordinate in zip(components, coordinates)
    )


def continuity_equation(
    density=rho_field, velocity=_DEFAULT_VELOCITY, time=t, coordinates=_COORDINATES
):
    r"""Return the compressible mass-conservation residual.

    A solution satisfies ``∂ρ/∂t + ∇·(ρ u) = 0``; this function returns
    the left-hand side, which should equal zero.
    """
    components = tuple(velocity)
    if len(components) != len(coordinates):
        raise ValueError("velocity and coordinate dimensions must match")
    mass_flux = tuple(density * component for component in components)
    return smp.diff(density, time) + divergence(mass_flux, coordinates)


def incompressible_continuity_equation(
    velocity=_DEFAULT_VELOCITY, coordinates=_COORDINATES
):
    r"""Return ``∇·u``, the incompressible continuity residual (equal to zero)."""
    return divergence(velocity, coordinates)


def navier_stokes_residual(
    density=rho,
    pressure=p,
    velocity=_DEFAULT_VELOCITY,
    dynamic_viscosity=mu,
    body_force=_DEFAULT_BODY_FORCE,
    time=t,
    coordinates=_COORDINATES,
):
    r"""Return the incompressible Newtonian momentum-equation residual.

    Uses constant density and constant dynamic viscosity:
    ``ρ(∂u/∂t + (u·∇)u) + ∇p - μ∇²u - ρf = 0``.
    Here ``f`` is body force per unit mass (acceleration), so gravity may be
    supplied as a component of ``body_force``. The function returns the
    three-component left-hand side; each component should equal zero.
    """
    velocity = tuple(velocity)
    body_force = tuple(body_force)
    if len(velocity) != len(coordinates) or len(body_force) != len(coordinates):
        raise ValueError(
            "velocity, body force, and coordinates must have matching dimensions"
        )

    acceleration = smp.Matrix(
        [
            material_derivative(component, velocity, time, coordinates)
            for component in velocity
        ]
    )
    pressure_gradient = gradient(pressure, coordinates)
    viscous_term = vector_laplacian(velocity, coordinates)
    return (
        density * acceleration
        + pressure_gradient
        - dynamic_viscosity * viscous_term
        - density * smp.Matrix(body_force)
    )


def reynolds_number(
    density=rho, speed=ux, characteristic_length=R, dynamic_viscosity=mu
):
    r"""Return ``Re = ρ U L / μ``, the dimensionless inertia-to-viscosity ratio."""
    return density * speed * characteristic_length / dynamic_viscosity


def kinematic_viscosity(dynamic_viscosity=mu, density=rho):
    r"""Return ``ν = μ/ρ`` (kinematic viscosity, with units of area per time)."""
    return dynamic_viscosity / density


def dynamic_pressure(density=rho, speed=ux):
    r"""Return ``q = ρ U²/2``, the incompressible dynamic pressure."""
    return density * speed**2 / 2


def hydrostatic_pressure(
    reference_pressure=p, density=rho, gravity=g, height=z, reference_height=0
):
    r"""Return pressure at height ``height`` in a static, constant-density fluid.

    The vertical coordinate is positive upward, giving ``p = p₀ - ρg(h-h₀)``.
    """
    return reference_pressure - density * gravity * (height - reference_height)


def volumetric_flow_rate(area=R**2 * smp.pi, average_speed=ux):
    r"""Return ``Q = A Ū`` for a cross-section with average normal speed ``Ū``."""
    return area * average_speed


def mass_flow_rate(density=rho, area=R**2 * smp.pi, average_speed=ux):
    r"""Return ``ṁ = ρ A Ū`` for uniform density and average normal speed."""
    return density * area * average_speed


def bernoulli_pressure(
    reference_pressure=p,
    density=rho,
    reference_speed=ux,
    speed=uy,
    gravity=g,
    reference_height=0,
    height=z,
):
    r"""Return pressure from Bernoulli's relation along a steady streamline.

    Assumes an incompressible, inviscid flow with no pump or turbine work:
    ``p + ρU²/2 + ρgh`` is constant along the streamline.
    """
    return (
        reference_pressure
        + density * (reference_speed**2 - speed**2) / 2
        + density * gravity * (reference_height - height)
    )


__all__ = [
    "velocity_vector",
    "gradient",
    "divergence",
    "laplacian",
    "vector_laplacian",
    "material_derivative",
    "continuity_equation",
    "incompressible_continuity_equation",
    "navier_stokes_residual",
    "reynolds_number",
    "kinematic_viscosity",
    "dynamic_pressure",
    "hydrostatic_pressure",
    "volumetric_flow_rate",
    "mass_flow_rate",
    "bernoulli_pressure",
]
