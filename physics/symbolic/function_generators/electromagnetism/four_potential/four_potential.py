import sympy as smp

from physics.symbolic.variables import (
    Bx,
    By,
    Bz,
    Ex,
    Ey,
    Ez,
    R,
    Sx,
    Sy,
    Sz,
    _path_parameter,
    c,
    eo,
    mo,
    omega,
    phi,
    q,
    t,
    x,
    y,
    z,
    vbar,
)

_OBSERVATION_POINT = (Sx, Sy, Sz)


def four_velocity(v=vbar, axis=1, dimension=4):
    """Return a velocity column in 3D or 4D."""
    if dimension not in (3, 4):
        raise ValueError("dimension must be 3 or 4")
    if not isinstance(axis, int) or not 1 <= axis < dimension:
        raise ValueError(f"axis must be between 1 and {dimension - 1}")

    velocity = [smp.S.One] + [smp.S.Zero] * (dimension - 1)
    velocity[axis] = v
    return smp.Matrix(velocity)


def _angular_integral(vector, angle):
    """Return a component-wise symbolic integral of a three-vector."""
    vector = smp.Matrix(vector)
    if vector.shape != (3, 1):
        raise ValueError("field integrands must be three-component column vectors")
    return smp.Matrix(
        [smp.Integral(component, (angle, 0, 2 * smp.pi)) for component in vector]
    )


def _ring_radiation_integrand(n_s, beta_s, beta_dot_s, kappa_s, R_s):
    """Return the retarded Liénard–Wiechert integrand for one ring element."""
    n_s = smp.Matrix(n_s)
    beta_s = smp.Matrix(beta_s)
    beta_dot_s = smp.Matrix(beta_dot_s)
    if any(vector.shape != (3, 1) for vector in (n_s, beta_s, beta_dot_s)):
        raise ValueError("n_s, beta_s, and beta_dot_s must be three-component columns")

    beta_squared = beta_s.dot(beta_s)
    velocity_term = (1 - beta_squared) * (n_s - beta_s) / (kappa_s**3 * R_s**2)
    acceleration_term = n_s.cross((n_s - beta_s).cross(beta_dot_s)) / (
        c * kappa_s**3 * R_s
    )
    return velocity_term + acceleration_term


def ring_electric_field(Q, n_s, beta_s, beta_dot_s, kappa_s, R_s, angle=None):
    r"""Return the symbolic electric field of a uniformly charged moving ring.

    Inputs describe each infinitesimal source element as the ring parameter
    ``angle`` varies. The caller must supply these quantities already evaluated
    at that element's retarded time ``t_r(angle)``:

    - ``Q``: total charge distributed uniformly around the ring, in coulombs.
    - ``n_s``: 3x1 unit vector pointing from the retarded source element to the
      observation point.
    - ``beta_s``: 3x1 dimensionless source velocity, ``v_s / c``.
    - ``beta_dot_s``: 3x1 time derivative of ``beta_s``, in inverse seconds.
    - ``kappa_s``: retarded factor ``1 - n_s.dot(beta_s)`` (often written
      ``kappa`` in Liénard–Wiechert formulas).
    - ``R_s``: source-to-observer distance at retarded time, in metres.
    - ``angle``: symbolic ring parameter integrated from 0 to ``2*pi``; by
      default this is a real SymPy symbol named ``varphi``.

    Returns a 3x1 SymPy matrix of component-wise unevaluated integrals. Pass
    each angle-dependent source quantity as a SymPy expression in ``angle``;
    this function does not solve the retarded-time equation itself.
    """
    if angle is None:
        angle = smp.Symbol("varphi", real=True)
    integrand = _ring_radiation_integrand(n_s, beta_s, beta_dot_s, kappa_s, R_s)
    return Q * _angular_integral(integrand, angle) / (8 * smp.pi**2 * eo)


def ring_magnetic_field(Q, n_s, beta_s, beta_dot_s, kappa_s, R_s, angle=None):
    r"""Return the symbolic magnetic field of a uniformly charged moving ring.

    Inputs have the same meanings and units as in
    :func:`ring_electric_field`: ``Q`` is total ring charge; ``n_s`` is the
    source-to-observer unit vector; ``beta_s = v_s/c``; ``beta_dot_s`` is its
    retarded-time derivative; ``kappa_s = 1 - n_s.dot(beta_s)``; ``R_s`` is
    retarded distance; and ``angle`` parameterizes the ring. Each input other
    than ``Q`` must be supplied with its retarded-time dependence on ``angle``.

    Returns a 3x1 SymPy matrix of component-wise unevaluated integrals of
    ``n_s`` crossed with the electric-field integrand, including the SI
    prefactor. Retarded time is not solved by this function.
    """
    if angle is None:
        angle = smp.Symbol("varphi", real=True)
    n_s = smp.Matrix(n_s)
    if n_s.shape != (3, 1):
        raise ValueError("n_s must be a three-component column vector")
    integrand = _ring_radiation_integrand(n_s, beta_s, beta_dot_s, kappa_s, R_s)
    return Q * _angular_integral(n_s.cross(integrand), angle) / (8 * smp.pi**2 * eo * c)


def _field_tensor(electric_field, magnetic_field):
    """Build the package's covariant field tensor from Cartesian E and B."""
    ex, ey, ez = electric_field
    bx, by, bz = magnetic_field
    return smp.Matrix(
        [
            [0, ex, ey, ez],
            [-ex, 0, bz, -by],
            [-ey, -bz, 0, bx],
            [-ez, by, -bx, 0],
        ]
    )


def _point_charge_fields(charge, source_position, velocity):
    """Return E and B for a uniformly moving point charge at an instant."""
    displacement = smp.Matrix(
        [point - source for point, source in zip(_OBSERVATION_POINT, source_position)]
    )
    velocity = smp.Matrix(velocity)
    beta = velocity / c
    beta_squared = beta.dot(beta)
    denominator = (
        displacement.dot(displacement)
        - beta.cross(displacement).dot(beta.cross(displacement))
    ) ** smp.Rational(3, 2)
    electric_field = (
        charge * (1 - beta_squared) * displacement / (4 * smp.pi * eo * denominator)
    )
    magnetic_field = velocity.cross(electric_field) / c**2
    return electric_field, magnetic_field


def _rotating_shell_fields(charge, radius, angular_velocity, position):
    """Return the rest-frame E and B fields of a uniformly charged shell."""
    displacement = smp.Matrix(position)
    angular_velocity = smp.Matrix(angular_velocity)
    distance_squared = displacement.dot(displacement)
    distance = smp.sqrt(distance_squared)
    magnetic_moment = charge * radius**2 * angular_velocity / 3

    electric_outside = charge * displacement / (4 * smp.pi * eo * distance**3)
    magnetic_outside = (
        mo
        * (
            3 * magnetic_moment.dot(displacement) * displacement
            - magnetic_moment * distance_squared
        )
        / (4 * smp.pi * distance**5)
    )
    electric_inside = smp.zeros(3, 1)
    magnetic_inside = mo * magnetic_moment / (2 * smp.pi * radius**3)

    electric_field = smp.Matrix(
        [
            smp.Piecewise(
                (electric_inside[i], distance < radius), (electric_outside[i], True)
            )
            for i in range(3)
        ]
    )
    magnetic_field = smp.Matrix(
        [
            smp.Piecewise(
                (magnetic_inside[i], distance < radius), (magnetic_outside[i], True)
            )
            for i in range(3)
        ]
    )
    return electric_field, magnetic_field


def four_potential(dimension=4):
    """Return the electromagnetic potential as a 3- or 4-component column."""
    if dimension not in (3, 4):
        raise ValueError("dimension must be 3 or 4")

    spatial_coordinates = (Sx, Sy, Sz)[: dimension - 1]
    component_names = ("A_x", "A_y", "A_z")[: dimension - 1]
    spatial_components = [
        smp.Function(name)(*spatial_coordinates) for name in component_names
    ]
    return smp.Matrix([phi / c, *spatial_components])


def electromagnetic_tensor(potential=None):
    """Return the field-strength tensor for a 3- or 4-component potential.

    The potential is ``(A_0, A_1, ...)``. Its spatial components depend on
    the available spatial coordinates. As in the original 4D function,
    temporal derivatives are omitted.
    """
    if potential is None:
        potential = four_potential()
    else:
        potential = smp.Matrix(potential)

    if potential.cols != 1 or potential.rows not in (3, 4):
        raise ValueError("potential must be a 3x1 or 4x1 column")

    dimension = potential.rows
    spatial_coordinates = (Sx, Sy, Sz)[: dimension - 1]
    tensor = smp.zeros(dimension, dimension)

    for index, coordinate in enumerate(spatial_coordinates, start=1):
        tensor[0, index] = -smp.diff(potential[0], coordinate)
        tensor[index, 0] = -tensor[0, index]

    for i, coordinate_i in enumerate(spatial_coordinates, start=1):
        for j, coordinate_j in enumerate(spatial_coordinates, start=1):
            tensor[i, j] = smp.diff(potential[j], coordinate_i) - smp.diff(
                potential[i], coordinate_j
            )

    return tensor


def four_potential_from_electromagnetic_tensor(
    electromagnetic_field_tensor,
    coordinates=None,
):
    r"""Recover a gauge-fixed potential from a closed 3x3 or 4x4 tensor.

    Uses Fock--Schwinger gauge,
    ``A_mu(x) = integral_0^1 s*x^nu*F_{nu mu}(s*x) ds``.
    Coordinates are ordered with time first; the time component in the path
    is ``ct``.
    """
    field_tensor = smp.Matrix(electromagnetic_field_tensor)
    if field_tensor.rows != field_tensor.cols or field_tensor.rows not in (3, 4):
        raise ValueError("electromagnetic_field_tensor must be 3x3 or 4x4")

    dimension = field_tensor.rows
    if coordinates is None:
        coordinates = (t, x, y, z)[:dimension]
    else:
        coordinates = tuple(coordinates)

    if len(coordinates) != dimension:
        raise ValueError("coordinates must match the tensor dimension")

    spacetime_coordinates = (c * coordinates[0], *coordinates[1:])
    path_substitutions = {
        coordinate: _path_parameter * coordinate for coordinate in coordinates
    }
    potential_components = [
        smp.integrate(
            _path_parameter
            * sum(
                spacetime_coordinates[nu]
                * field_tensor[nu, mu].subs(path_substitutions, simultaneous=True)
                for nu in range(dimension)
            ),
            (_path_parameter, 0, 1),
        )
        for mu in range(dimension)
    ]
    return smp.Matrix(potential_components)


def lorentz_four_force(
    electromagnetic_tensor,
    four_velocity,
    charge=q,
    index_position="upper",
    metric=None,
):
    r"""Return the contravariant Lorentz force in 3 or 4 dimensions.

    The field tensor, velocity, and metric must have matching dimensions.
    The default metric has signature ``(-,+,+)`` in 3D and
    ``(-,+,+,+)`` in 4D.
    """
    field_tensor = smp.Matrix(electromagnetic_tensor)
    velocity = smp.Matrix(four_velocity)

    if field_tensor.rows != field_tensor.cols or field_tensor.rows not in (3, 4):
        raise ValueError("electromagnetic_tensor must be 3x3 or 4x4")

    dimension = field_tensor.rows
    if velocity.shape == (1, dimension):
        velocity = velocity.T
    if velocity.shape != (dimension, 1):
        raise ValueError(f"four_velocity must contain {dimension} components")

    if metric is None:
        metric = smp.diag(-1, *([1] * (dimension - 1)))
    else:
        metric = smp.Matrix(metric)

    if metric.shape != (dimension, dimension):
        raise ValueError(f"metric must be a {dimension}x{dimension} matrix")

    if index_position in ("upper", "superscript", "super"):
        return charge * field_tensor * metric * velocity
    if index_position in ("lower", "subscript", "sub"):
        return charge * metric.inv() * field_tensor * velocity
    raise ValueError(
        "index_position must be 'upper'/'superscript' or 'lower'/'subscript'"
    )


def stationary_point_charge_tensor(charge=None, source_position=(0, 0, 0)):
    """Return F for a stationary point charge in SI units.

    The observation point is ``(Sx, Sy, Sz)`` and ``source_position`` is the
    charge location. The tensor uses the same component convention as
    :func:`physics.symbolic.tensors_kinda.F.F_sub_mu_nu`.
    """
    if charge is None:
        charge = q
    electric_field, magnetic_field = _point_charge_fields(
        charge, source_position, (0, 0, 0)
    )
    return _field_tensor(electric_field, magnetic_field)


def moving_point_charge_tensor(charge=None, velocity=None, source_position=(0, 0, 0)):
    """Return the constant-velocity point-charge field tensor in SI units.

    This is the uniform-motion (velocity-field) solution evaluated at an
    instant, with displacement from ``source_position`` to ``(Sx, Sy, Sz)``.
    Radiation and retardation from acceleration are not included.
    """
    if charge is None:
        charge = q
    if velocity is None:
        velocity = tuple(smp.Symbol(f"v_{axis}") for axis in "xyz")
    electric_field, magnetic_field = _point_charge_fields(
        charge, source_position, velocity
    )
    return _field_tensor(electric_field, magnetic_field)


def co_moving_charges_tensor(charges, velocity=None, source_positions=None):
    """Return F for point charges sharing one constant velocity.

    ``charges`` and ``source_positions`` are sequences of charges and their
    instantaneous positions. Fields are summed by linear superposition.
    """
    charges = tuple(charges)
    if source_positions is None:
        source_positions = ((0, 0, 0),) * len(charges)
    source_positions = tuple(source_positions)
    if len(charges) != len(source_positions):
        raise ValueError("charges and source_positions must have equal lengths")
    if velocity is None:
        velocity = tuple(smp.Symbol(f"v_{axis}") for axis in "xyz")

    electric_field = smp.zeros(3, 1)
    magnetic_field = smp.zeros(3, 1)
    for charge, position in zip(charges, source_positions):
        charge_electric, charge_magnetic = _point_charge_fields(
            charge, position, velocity
        )
        electric_field += charge_electric
        magnetic_field += charge_magnetic
    return _field_tensor(electric_field, magnetic_field)


def stationary_rotating_spherical_charge_tensor(
    charge=None, radius=None, angular_velocity=None
):
    """Return F for a uniformly charged, stationary, rotating thin shell.

    The shell has uniform surface charge, radius ``radius``, and angular
    velocity vector ``angular_velocity``. Its magnetic moment is
    ``charge * radius**2 * angular_velocity / 3``. The observation point is
    ``(Sx, Sy, Sz)`` relative to the shell center. The fields are piecewise:
    Coulomb and dipole fields outside, zero electric field and uniform
    magnetic field inside. The idealized field is discontinuous at the shell.
    """
    if charge is None:
        charge = q
    if radius is None:
        radius = R
    if angular_velocity is None:
        angular_velocity = tuple(smp.Symbol(f"omega_{axis}") for axis in "xyz")
    electric_field, magnetic_field = _rotating_shell_fields(
        charge, radius, angular_velocity, _OBSERVATION_POINT
    )
    return _field_tensor(electric_field, magnetic_field)


def moving_rotating_spherical_charge_tensor(
    charge=None, radius=None, angular_velocity=None, velocity=None
):
    """Return F for a uniformly rotating charged shell in uniform translation.

    The shell model matches :func:`stationary_rotating_spherical_charge_tensor`.
    ``velocity`` is the shell's translational velocity in the lab frame. The
    rest-frame static fields are evaluated at Lorentz-transformed coordinates
    and the electric and magnetic fields are Lorentz-boosted into the lab
    frame. Acceleration and radiation are not included.
    """
    if charge is None:
        charge = q
    if radius is None:
        radius = R
    if angular_velocity is None:
        angular_velocity = tuple(smp.Symbol(f"omega_{axis}") for axis in "xyz")
    if velocity is None:
        velocity = tuple(smp.Symbol(f"u_{axis}") for axis in "xyz")

    velocity = smp.Matrix(velocity)
    velocity_squared = velocity.dot(velocity)
    gamma = 1 / smp.sqrt(1 - velocity_squared / c**2)
    boost_factor = gamma**2 / ((gamma + 1) * c**2)
    lab_position = smp.Matrix(_OBSERVATION_POINT)
    rest_position = (
        lab_position
        + (boost_factor * velocity.dot(lab_position) - gamma * t) * velocity
    )

    rest_electric, rest_magnetic = _rotating_shell_fields(
        charge, radius, angular_velocity, rest_position
    )
    electric_field = (
        gamma * (rest_electric - velocity.cross(rest_magnetic))
        - boost_factor * velocity.dot(rest_electric) * velocity
    )
    magnetic_field = (
        gamma * (rest_magnetic + velocity.cross(rest_electric) / c**2)
        - boost_factor * velocity.dot(rest_magnetic) * velocity
    )
    return _field_tensor(electric_field, magnetic_field)


if __name__ == "__main__":
    # Symbolic uniformly rotating ring: its charge, radius, angular speed,
    # observation coordinates, and light speed come from variables.py.
    ring_angle = smp.Symbol("ring_angle", real=True)
    ring_speed_fraction = omega * R / c
    source_velocity = smp.Matrix(
        [
            -ring_speed_fraction * smp.sin(ring_angle),
            ring_speed_fraction * smp.cos(ring_angle),
            0,
        ]
    )
    source_acceleration = smp.Matrix(
        [
            -omega * ring_speed_fraction * smp.cos(ring_angle),
            -omega * ring_speed_fraction * smp.sin(ring_angle),
            0,
        ]
    )
    displacement = smp.Matrix(
        [
            Sx - R * smp.cos(ring_angle),
            Sy - R * smp.sin(ring_angle),
            Sz,
        ]
    )
    source_to_observer_distance = smp.sqrt(displacement.dot(displacement))
    source_to_observer_direction = displacement / source_to_observer_distance
    retarded_factor = 1 - source_to_observer_direction.dot(source_velocity)
    ring_electric_field(
        q,
        source_to_observer_direction,
        source_velocity,
        source_acceleration,
        retarded_factor,
        source_to_observer_distance,
        ring_angle,
    )
    ring_magnetic_field(
        q,
        source_to_observer_direction,
        source_velocity,
        source_acceleration,
        retarded_factor,
        source_to_observer_distance,
        ring_angle,
    )
