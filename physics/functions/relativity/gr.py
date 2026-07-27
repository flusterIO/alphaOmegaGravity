import numpy as np
import scipy.constants as spc
import astropy.constants.iau2015 as ac


def schwartzchildRadius(mass, G=spc.G, c=spc.c):
    return 2 * G * mass / c**2


def schwartzchildMetric(mass=ac.M_earth.value, G=spc.G, r=ac.R_earth.value):
    return 1 - 2 * G * mass / r


def GR_satellite_radius(r_sat):
    r"""
    Computes the scalar Euclidean distance of the satellite from the Earth's center of mass.

    LaTeX:
    $$
    r_s = \sqrt{x_s^2 + y_s^2 + z_s^2}
    $$

    Description:
    This represents the instantaneous orbital radius in the Earth-Centered Inertial (ECI)
    coordinate frame. It is required to determine the Earth's local gravitational
    potential at the satellite's exact position in space.

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite Cartesian position vector [x, y, z] in meters (ECI frame).

    Returns:
    --------
    np.float64 or ndarray
        Scalar satellite radius in meters.
    """
    r_sat = np.asarray(r_sat, dtype=np.float64)
    return np.linalg.norm(r_sat, axis=-1)


def GR_satellite_velocity_squared(v_sat):
    r"""
    Computes the scalar velocity squared of the satellite in the non-rotating ECI frame.

    LaTeX:
    $$
    v_s^2 = v_{xs}^2 + v_{ys}^2 + v_{zs}^2
    $$

    Description:
    This square of the orbital speed determines the second-order Doppler time dilation
    (kinematic redshift) experienced by the moving satellite clock relative to a
    hypothetical stationary clock at infinity.

    Parameters:
    -----------
    v_sat : array_like, shape (3,) or (N, 3)
        Satellite Cartesian velocity vector [vx, vy, vz] in m/s (ECI frame).

    Returns:
    --------
    np.float64 or ndarray
        Scalar velocity squared in (m/s)^2.
    """
    v_sat = np.asarray(v_sat, dtype=np.float64)
    return np.sum(v_sat**2, axis=-1)


def GR_receiver_radius(r_rec):
    r"""
    Computes the scalar Euclidean distance of the ground receiver from the Earth's center.

    LaTeX:
    $$
    r_r = \sqrt{x_r^2 + y_r^2 + z_r^2}
    $$

    Description:
    This defines the geocentric distance to the user's receiving antenna. It is used in
    conjunction with the satellite radius to evaluate line-of-sight gravitational
    retardation (Shapiro delay).

    Parameters:
    -----------
    r_rec : array_like, shape (3,) or (N, 3)
        Receiver Cartesian position vector [x, y, z] in meters (ECI frame).

    Returns:
    --------
    np.float64 or ndarray
        Scalar receiver radius in meters.
    """
    r_rec = np.asarray(r_rec, dtype=np.float64)
    return np.linalg.norm(r_rec, axis=-1)


def GR_geometric_distance(r_sat, r_rec):
    r"""
    Computes the straight-line Euclidean distance between transmitter and receiver.

    LaTeX:
    $$
    R_{sr} = \sqrt{(x_s - x_r)^2 + (y_s - y_r)^2 + (z_s - z_r)^2}
    $$

    Description:
    This is the uncorrected geometric range that the radio photon traverses through vacuum.
    Dividing this by the speed of light gives the baseline coordinate propagation time,
    to which relativistic path delays are added.

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite position vector [x, y, z] in meters (ECI frame).
    r_rec : array_like, shape (3,) or (N, 3)
        Receiver position vector [x, y, z] in meters (ECI frame).

    Returns:
    --------
    np.float64 or ndarray
        Geometric range R_sr in meters.
    """
    r_sat = np.asarray(r_sat, dtype=np.float64)
    r_rec = np.asarray(r_rec, dtype=np.float64)
    return np.linalg.norm(r_sat - r_rec, axis=-1)


def GR_semi_major_axis(r_sat, v_sat, GM=3.986004418e14):
    r"""
    Derives the orbital semi-major axis dynamically from Cartesian state vectors.

    LaTeX:
    $$
    a = \left( \frac{2}{r_s} - \frac{v_s^2}{GM} \right)^{-1}
    $$

    Description:
    Using the conservation of specific orbital energy (vis-viva equation), this extracts
    the constant semi-major axis of the orbit directly from instantaneous position and
    velocity without requiring external Keplerian ephemeris data.

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite position vector [x, y, z] in meters (ECI frame).
    v_sat : array_like, shape (3,) or (N, 3)
        Satellite velocity vector [vx, vy, vz] in m/s (ECI frame).
    GM : float, optional
        Earth's gravitational parameter in m^3/s^2 (default: WGS84 value).

    Returns:
    --------
    np.float64 or ndarray
        Semi-major axis 'a' in meters.
    """
    r_s = GR_satellite_radius(r_sat)
    v_s_sq = GR_satellite_velocity_squared(v_sat)
    return 1.0 / (2.0 / r_s - v_s_sq / GM)


def GR_secular_fractional_offset(
    a, GM=3.986004418e14, c=299792458.0, Phi_0=-62636856.0
):
    r"""
    Computes the steady-state dimensionless frequency offset of the satellite clock.

    LaTeX:
    $$
    \frac{\Delta f}{f_0} = \frac{3GM}{2ac^2} + \frac{\Phi_0}{c^2}
    $$

    Description:
    This combines the gravitational blueshift (higher altitude clock ticking faster) and
    the kinematic redshift (moving clock ticking slower) into a single net fractional
    frequency difference relative to a reference clock resting on Earth's geoid (Phi_0).
    For nominal GPS orbits, this resolves to approximately +4.4647e-10.

    Parameters:
    -----------
    a : float or array_like
        Orbital semi-major axis in meters.
    GM : float, optional
        Earth's gravitational parameter in m^3/s^2.
    c : float, optional
        Speed of light in vacuum in m/s.
    Phi_0 : float, optional
        Reference gravitational potential on Earth's geoid in m^2/s^2.

    Returns:
    --------
    np.float64 or ndarray
        Dimensionless fractional frequency offset (df / f0).
    """
    a = np.asarray(a, dtype=np.float64)
    return (3.0 * GM) / (2.0 * a * (c**2)) + (Phi_0 / (c**2))


def GR_daily_secular_drift(secular_offset):
    r"""
    Converts a fractional frequency offset into daily accumulated time drift in microseconds.

    LaTeX:
    $$
    \Delta t_{\text{secular}} = \left( \frac{\Delta f}{f_0} \right) \times 86400 \times 10^6
    $$

    Description:
    Translates the dimensionless rate difference into a tangible time drift over a standard
    24-hour solar day (86,400 seconds). For GPS, this yields the factory pre-shift target
    of +38.58 microseconds per day that is baked into satellite synthesizers before launch.

    Parameters:
    -----------
    secular_offset : float or array_like
        Dimensionless fractional frequency offset (df / f0).

    Returns:
    --------
    np.float64 or ndarray
        Accumulated time drift in microseconds per day (us/day).
    """
    secular_offset = np.asarray(secular_offset, dtype=np.float64)
    return secular_offset * 86400.0 * 1e6


def GR_eccentricity_correction(r_sat, v_sat, c=299792458.0):
    r"""
    Computes the periodic proper time correction due to orbital eccentricity.

    LaTeX:
    $$
    \Delta t_{\text{ecc}} = -\frac{2}{c^2} (x_s v_{xs} + y_s v_{ys} + z_s v_{zs})
    $$

    Description:
    Because elliptical orbits cause continuous velocity and altitude oscillations, this
    dot-product term isolates and corrects the periodic clock error as the satellite
    transits between perigee and apogee. This correction is evaluated continuously by
    the GPS receiver software.

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite position vector [x, y, z] in meters (ECI frame).
    v_sat : array_like, shape (3,) or (N, 3)
        Satellite velocity vector [vx, vy, vz] in m/s (ECI frame).
    c : float, optional
        Speed of light in vacuum in m/s.

    Returns:
    --------
    np.float64 or ndarray
        Eccentricity time correction in seconds.
    """
    r_sat = np.asarray(r_sat, dtype=np.float64)
    v_sat = np.asarray(v_sat, dtype=np.float64)
    r_dot_v = np.sum(r_sat * v_sat, axis=-1)
    return -(2.0 / (c**2)) * r_dot_v


def GR_sagnac_delay(r_sat, r_rec, omega_e=7.2921151467e-5, c=299792458.0):
    r"""
    Computes the Sagnac time delay accounting for Earth rotation during signal transit.

    LaTeX:
    $$
    \Delta t_{\text{Sagnac}} = \frac{\omega_e (x_s y_r - y_s x_r)}{c^2}
    $$

    Description:
    During the ~70 milliseconds a photon travels from satellite to Earth, an Earth-fixed
    receiver rotates eastward away from or toward the signal path. This simplified
    cross-product expansion (valid because Earth rotates solely around the Z-axis) corrects
    the ECI propagation time to the rotating Earth-Centered, Earth-Fixed (ECEF) user frame.

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite position vector [x, y, z] in meters (ECI frame).
    r_rec : array_like, shape (3,) or (N, 3)
        Receiver position vector [x, y, z] in meters (ECI frame).
    omega_e : float, optional
        Earth angular rotation rate in rad/s (default: WGS84 value).
    c : float, optional
        Speed of light in vacuum in m/s.

    Returns:
    --------
    np.float64 or ndarray
        Sagnac signal transit time correction in seconds.
    """
    r_sat = np.asarray(r_sat, dtype=np.float64)
    r_rec = np.asarray(r_rec, dtype=np.float64)

    if r_sat.ndim == 1:
        sagnac_num = r_sat[0] * r_rec[1] - r_sat[1] * r_rec[0]
    else:
        sagnac_num = r_sat[:, 0] * r_rec[:, 1] - r_sat[:, 1] * r_rec[:, 0]

    return (omega_e * sagnac_num) / (c**2)


def GR_shapiro_delay(r_sat, r_rec, GM=3.986004418e14, c=299792458.0):
    r"""
    Computes the Shapiro line-of-sight gravitational signal retardation.

    LaTeX:
    $$
    \Delta t_{\text{Shapiro}} = \frac{2GM}{c^3} \ln \left( \frac{r_s + r_r + R_{sr}}{r_s + r_r - R_{sr}} \right)
    $$

    Description:
    According to general relativity, electromagnetic signals dip into Earth's gravitational
    potential well along their flight path, experiencing a slight coordinate time delay.
    This logarithmic integration along the scalar path length adds up to ~180 picoseconds
    of delay (equivalent to roughly 2 centimeters of pseudorange error).

    Parameters:
    -----------
    r_sat : array_like, shape (3,) or (N, 3)
        Satellite position vector [x, y, z] in meters (ECI frame).
    r_rec : array_like, shape (3,) or (N, 3)
        Receiver position vector [x, y, z] in meters (ECI frame).
    GM : float, optional
        Earth's gravitational parameter in m^3/s^2.
    c : float, optional
        Speed of light in vacuum in m/s.

    Returns:
    --------
    np.float64 or ndarray
        Shapiro gravitational time delay in seconds.
    """
    r_s = GR_satellite_radius(r_sat)
    r_r = GR_receiver_radius(r_rec)
    R_sr = GR_geometric_distance(r_sat, r_rec)

    shapiro_arg = (r_s + r_r + R_sr) / (r_s + r_r - R_sr)
    return (2.0 * GM / (c**3)) * np.log(shapiro_arg)
