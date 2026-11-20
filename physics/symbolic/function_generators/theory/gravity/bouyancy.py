"""Symbolic equations for buoyancy and hydrostatic pressure."""

import sympy as smp

from physics.symbolic.variables import g, m, rho, depth, V, rho_object


def buoyant_force(fluid_density=rho, gravity=g, displaced_volume=V):
    """Return Archimedes' buoyant force, positive upward: rho_fluid * g * V."""
    return fluid_density * gravity * displaced_volume


def object_weight(mass=m, gravity=g):
    """Return an object's weight magnitude: m * g."""
    return mass * gravity


def net_buoyant_force(
    fluid_density=rho,
    gravity=g,
    displaced_volume=V,
    object_mass=m,
):
    """Return the net vertical force, positive upward."""
    return buoyant_force(fluid_density, gravity, displaced_volume) - object_weight(
        object_mass, gravity
    )


def apparent_weight(
    object_mass=m,
    gravity=g,
    fluid_density=rho,
    displaced_volume=V,
):
    """Return the object's apparent weight while immersed in a fluid."""
    return object_weight(object_mass, gravity) - buoyant_force(
        fluid_density, gravity, displaced_volume
    )


def hydrostatic_pressure(fluid_density=rho, gravity=g, fluid_depth=depth):
    """Return gauge pressure at depth: rho_fluid * g * depth."""
    return fluid_density * gravity * fluid_depth


def displaced_volume_for_floating_object(object_mass=m, fluid_density=rho):
    """Return the displaced volume required for floating equilibrium: m/rho."""
    return object_mass / fluid_density


def submerged_fraction(object_density=rho_object, fluid_density=rho):
    """Return the submerged fraction of a floating uniform object."""
    return object_density / fluid_density


def floats_or_sinks_condition(object_density=rho_object, fluid_density=rho):
    """Return the density difference; negative floats, positive sinks."""
    return object_density - fluid_density
