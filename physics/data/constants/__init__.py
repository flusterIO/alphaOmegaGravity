"""Supplementary cosmology and large-number estimates.

Standard physical and astronomical constants are imported from SciPy and
Astropy rather than redefined here. Public names in this module are derived
scales, counts, or observational/model-dependent estimates. Values are floats
in SI units unless a name says otherwise.

Requires SciPy and Astropy. Cosmology uses Astropy's Planck18 model.

Earth atmospheric-electricity estimates are from Harrison et al. (2020),
"Extensive layer clouds in the global electric circuit: their effects on
vertical charge distribution and storage," Proceedings of the Royal Society A
476, 20190758, https://doi.org/10.1098/rspa.2019.0758. Their typical fair-weather
surface field of -120 V/m implies a surface charge density of about -1 nC/m².
Extrapolating that density over the whole Earth gives an estimated fair-weather
surface charge of -5×10⁵ C; the paper explicitly cautions that this is not the
total charge on Earth's surface. The paper also gives a typical upper (often
called ionospheric) potential of about 240 kV relative to the surface.
"""

from math import pi as _pi

from astropy import constants as _ac
from astropy.cosmology import Planck18 as _cosmo
from scipy import constants as _sc
from scipy.special import zeta as _zeta
import physics.data.li_observation as li
import physics.data.constants.milky_way as mw
import physics.data.constants.earth as earth
import physics.data.constants.GR as GR
import physics.data.eilers_observation as eilers

# Atmospheric electric quantities for Earth's fair-weather global circuit.
# These are estimates from Harrison et al. (2020), cited in this module's
# docstring; the integrated charge is not Earth's total surface charge.
EARTH_FAIR_WEATHER_ELECTRIC_FIELD_V_M = -120.0
EARTH_FAIR_WEATHER_SURFACE_CHARGE_DENSITY_C_M2 = -1.0e-9
EARTH_FAIR_WEATHER_SURFACE_CHARGE_ESTIMATE_C = -5.0e5
EARTH_IONOSPHERIC_POTENTIAL_V = 2.4e5

# Cosmological distances and time scales from Astropy's Planck 2018 model.
HUBBLE_PARAMETER_TODAY = _cosmo.H0.to_value("1 / s")
HUBBLE_TIME = 1.0 / HUBBLE_PARAMETER_TODAY
HUBBLE_LENGTH = (_ac.c / _cosmo.H0).to_value("m")  # c/H0, not particle-horizon radius
UNIVERSE_AGE = _cosmo.age(0).to_value("s")
UNIVERSE_AGE_GYR = _cosmo.age(0).to_value("Gyr")
# Comoving distance to very high redshift approximates the particle horizon.
OBSERVABLE_UNIVERSE_RADIUS = _cosmo.comoving_distance(1.0e5).to_value("m")
OBSERVABLE_UNIVERSE_DIAMETER = 2.0 * OBSERVABLE_UNIVERSE_RADIUS
OBSERVABLE_UNIVERSE_VOLUME = (4.0 / 3.0) * _pi * OBSERVABLE_UNIVERSE_RADIUS**3

# Present CMB photon density from the blackbody number-density formula
# n_gamma = 2*zeta(3)/pi^2 * (kT/(hbar*c))^3.
_T_CMB = _cosmo.Tcmb0.to_value("K")
PHOTON_DENSITY_CMB = (
    2.0 * _zeta(3.0) / _pi**2 * (_sc.Boltzmann * _T_CMB / (_sc.hbar * _sc.c)) ** 3
)
# Standard thermal relic neutrino-to-photon number ratio: 9/11 for three
# active species, including antiparticles, after electron-positron annihilation.
NEUTRINO_DENSITY_RELIC = (9.0 / 11.0) * PHOTON_DENSITY_CMB

# Present-day mass in the observable volume from Planck18 density fractions.
_RHO_CRITICAL = _cosmo.critical_density0.to_value("kg / m3")
_RHO_BARYON = _cosmo.Ob0 * _RHO_CRITICAL
_RHO_DARK_MATTER = (_cosmo.Om0 - _cosmo.Ob0) * _RHO_CRITICAL
_RHO_DARK_ENERGY = _cosmo.Ode0 * _RHO_CRITICAL

BARYON_MASS_OBSERVABLE = _RHO_BARYON * OBSERVABLE_UNIVERSE_VOLUME
DARK_MATTER_MASS_OBSERVABLE = _RHO_DARK_MATTER * OBSERVABLE_UNIVERSE_VOLUME
TOTAL_MATTER_MASS_OBSERVABLE = _cosmo.Om0 * _RHO_CRITICAL * OBSERVABLE_UNIVERSE_VOLUME
DARK_ENERGY_MASS_EQUIVALENT_OBSERVABLE = _RHO_DARK_ENERGY * OBSERVABLE_UNIVERSE_VOLUME
TOTAL_MATTER_ENERGY_OBSERVABLE = TOTAL_MATTER_MASS_OBSERVABLE * _ac.c.value**2

# Number counts use these model-derived masses, particle masses, and densities.
BARYON_COUNT_OBSERVABLE = BARYON_MASS_OBSERVABLE / _ac.m_p.value
PHOTON_COUNT_OBSERVABLE = PHOTON_DENSITY_CMB * OBSERVABLE_UNIVERSE_VOLUME
NEUTRINO_COUNT_OBSERVABLE = NEUTRINO_DENSITY_RELIC * OBSERVABLE_UNIVERSE_VOLUME
PHOTONS_PER_BARYON = PHOTON_COUNT_OBSERVABLE / BARYON_COUNT_OBSERVABLE
EDDINGTON_BARYON_NUMBER_SCALE = BARYON_COUNT_OBSERVABLE

# Characteristic scales derived from package constants.
PROTON_REST_ENERGY = (_ac.m_p * _ac.c**2).to_value("J")
PROTON_REST_ENERGY_EV = (_ac.m_p * _ac.c**2).to_value("eV")
PLANCK_LENGTH = (_ac.hbar * _ac.G / _ac.c**3) ** 0.5
PLANCK_LENGTH = PLANCK_LENGTH.to_value("m")
PLANCK_MASS = (_ac.hbar * _ac.c / _ac.G) ** 0.5
PLANCK_MASS = PLANCK_MASS.to_value("kg")
PLANCK_TIME = PLANCK_LENGTH / _ac.c.to_value("m / s")
ELECTRON_COMPTON_WAVELENGTH = (_ac.h / (_ac.m_e * _ac.c)).to_value("m")

# Gravity/electromagnetic force ratio for a proton-electron pair.
# Coulomb's constant is 1/(4*pi*epsilon_0), from the SciPy SI constants.
_COULOMB_K = 1.0 / (4.0 * _pi * _sc.epsilon_0)
PROTON_ELECTRON_GRAVITY_TO_ELECTROMAGNETIC_FORCE_RATIO = (
    _sc.G * _ac.m_p.value * _ac.m_e.value / (_COULOMB_K * _sc.e**2)
)
PROTON_ELECTRON_ELECTROMAGNETIC_TO_GRAVITY_FORCE_RATIO = (
    1.0 / PROTON_ELECTRON_GRAVITY_TO_ELECTROMAGNETIC_FORCE_RATIO
)
UNIVERSE_AGE_IN_PLANCK_TIMES = UNIVERSE_AGE / PLANCK_TIME

# Observed dark-energy density compared with the Planck density. The inverse
# is a data-derived scale related to (but not identical to) the vacuum puzzle.
_RHO_PLANCK = _ac.m_p.value / PLANCK_LENGTH**3
_RHO_PLANCK = PLANCK_MASS / PLANCK_LENGTH**3
DARK_ENERGY_TO_PLANCK_DENSITY_RATIO = _RHO_DARK_ENERGY / _RHO_PLANCK
PLANCK_TO_DARK_ENERGY_DENSITY_RATIO = 1.0 / DARK_ENERGY_TO_PLANCK_DENSITY_RATIO
