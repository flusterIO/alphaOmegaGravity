"""Unit conversion helpers for numeric physics calculations."""

import astropy.units as u


def gev_to_kg(energy_gev):
    """Convert energy in gigaelectronvolts to its mass equivalent in kilograms.

    ``energy_gev`` may be a scalar or array-like value. The returned value is
    in kilograms and uses Astropy's mass-energy equivalency and unit constants.
    """
    return u.Quantity(energy_gev, u.GeV).to_value(
        u.kg, equivalencies=u.mass_energy()
    )
