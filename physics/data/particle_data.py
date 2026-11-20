"""Build an extensible DataFrame of Standard Model elementary particles.

Masses and electric charges come from the 2026 Particle Data Group summary
tables:

- Quarks: https://pdg.lbl.gov/2026/tables/rpp2026-sum-quarks.pdf
- Leptons: https://pdg.lbl.gov/2026/tables/rpp2026-sum-leptons.pdf
- Gauge and Higgs bosons: https://pdg.lbl.gov/2026/tables/rpp2026-sum-gauge-higgs-bosons.pdf

Masses are expressed in GeV/c² and electric charges in elementary-charge
units and coulombs. Quark masses use the scale and scheme specified by PDG;
the top mass is the direct event-kinematics value. Neutrino flavor masses
remain undetermined, so their mass values are missing rather than set to zero.
The photon mass is also stored as a limit, not as a measured nonzero mass.
"""

from collections.abc import Mapping

import pandas as pd
from scipy.constants import c, elementary_charge, electron_volt


_GEV_C2_TO_KG = 1e9 * electron_volt / c**2


_PDG_QUARKS = "https://pdg.lbl.gov/2026/tables/rpp2026-sum-quarks.pdf"
_PDG_LEPTONS = "https://pdg.lbl.gov/2026/tables/rpp2026-sum-leptons.pdf"
_PDG_BOSONS = "https://pdg.lbl.gov/2026/tables/rpp2026-sum-gauge-higgs-bosons.pdf"


def _particle(
    symbol,
    name,
    category,
    charge_e,
    mass_gev_c2=None,
    mass_uncertainty_gev_c2=None,
    mass_upper_limit_gev_c2=None,
    mass_status="measured",
    mass_definition=None,
    source=None,
    antiparticle=False,
):
    """Create a normalized, extensible particle record."""
    return {
        "symbol": symbol,
        "name": name,
        "category": category,
        "is_antiparticle": antiparticle,
        "electric_charge_e": charge_e,
        "electric_charge_c": charge_e * elementary_charge,
        "mass_gev_c2": mass_gev_c2,
        "mass_kg": (mass_gev_c2 * _GEV_C2_TO_KG if mass_gev_c2 is not None else None),
        "mass_uncertainty_gev_c2": mass_uncertainty_gev_c2,
        "mass_upper_limit_gev_c2": mass_upper_limit_gev_c2,
        "mass_status": mass_status,
        "mass_definition": mass_definition,
        "source": source,
    }


def particleDataFrame(additional_data=None, source=False):
    """Return a DataFrame of Standard Model elementary-particle properties.

    The table includes the six quark flavors and their antiparticles, charged
    leptons and antileptons, neutrino and antineutrino flavors, and the photon,
    eight gluon states (represented by one row with ``state_count=8``), both W
    charge states, the Z boson, and the Higgs boson.

    ``additional_data`` may be a mapping from particle ``symbol`` to mappings
    of extra column names and values. For example,
    ``{"e-": {"generation": 1, "spin": 0.5}}`` adds those columns and
    populates the electron row, leaving other rows missing. Unknown symbols
    raise ``ValueError`` to avoid silently attaching data to the wrong row.

    Masses and charges are based on the 2026 Particle Data Group Summary
    Tables: quarks
    (https://pdg.lbl.gov/2026/tables/rpp2026-sum-quarks.pdf), leptons
    (https://pdg.lbl.gov/2026/tables/rpp2026-sum-leptons.pdf), and gauge/Higgs
    bosons
    (https://pdg.lbl.gov/2026/tables/rpp2026-sum-gauge-higgs-bosons.pdf).
    Masses are in GeV/c² and kilograms (``mass_kg`` is calculated from the
    nominal ``mass_gev_c2`` value); uncertainties and upper limits use GeV/c².
    Charges are provided both in units of the elementary charge and coulombs.
    """
    rows = []

    # PDG quark masses; u, d, s are MS at 2 GeV, c and b are MS at their
    # respective mass scales, and t is the direct event-kinematics mass.
    quarks = (
        ("u", "up", 2 / 3, 0.00216, 0.00007, "MS at 2 GeV"),
        ("d", "down", -1 / 3, 0.00470, 0.00007, "MS at 2 GeV"),
        ("s", "strange", -1 / 3, 0.0929, 0.0007, "MS at 2 GeV"),
        ("c", "charm", 2 / 3, 1.2729, 0.0045, "MS at the quark mass"),
        ("b", "bottom", -1 / 3, 4.186, 0.006, "MS at the quark mass"),
        ("t", "top", 2 / 3, 172.60, 0.27, "direct event-kinematics mass"),
    )
    for symbol, name, charge, mass, uncertainty, definition in quarks:
        rows.append(
            _particle(
                symbol,
                f"{name} quark",
                "quark",
                charge,
                mass,
                uncertainty,
                mass_definition=definition,
                source=_PDG_QUARKS,
            )
        )
        rows.append(
            _particle(
                f"anti_{symbol}",
                f"anti-{name} quark",
                "quark",
                -charge,
                mass,
                uncertainty,
                mass_definition=definition,
                source=_PDG_QUARKS,
                antiparticle=True,
            )
        )

    # PDG charged-lepton values, converted from MeV/c² to GeV/c².
    charged_leptons = (
        ("e-", "electron", -1, 0.51099895069e-3, 0.00000000016e-3),
        ("mu-", "muon", -1, 105.6583755e-3, 0.0000023e-3),
        ("tau-", "tau", -1, 1.77693, 0.00009),
    )
    for symbol, name, charge, mass, uncertainty in charged_leptons:
        rows.append(
            _particle(
                symbol,
                name,
                "charged lepton",
                charge,
                mass,
                uncertainty,
                source=_PDG_LEPTONS,
            )
        )
        antiparticle_symbol = symbol.replace("-", "+")
        rows.append(
            _particle(
                antiparticle_symbol,
                f"anti-{name}",
                "charged lepton",
                -charge,
                mass,
                uncertainty,
                source=_PDG_LEPTONS,
                antiparticle=True,
            )
        )

    # PDG gives a 90% CL direct tritium-decay bound on effective electron-
    # neutrino mass; flavor masses are not individually measured.
    electron_neutrino_limit = 0.8e-9  # 0.8 eV/c² expressed in GeV/c².
    for flavor in ("e", "mu", "tau"):
        for antiparticle, prefix in ((False, "nu_"), (True, "anti_nu_")):
            symbol = f"{prefix}{flavor}"
            is_electron_flavor = flavor == "e"
            rows.append(
                _particle(
                    symbol,
                    f"{'anti-' if antiparticle else ''}{flavor} neutrino",
                    "neutrino",
                    0,
                    mass_upper_limit_gev_c2=(
                        electron_neutrino_limit if is_electron_flavor else None
                    ),
                    mass_status=(
                        "effective_mass_upper_limit"
                        if is_electron_flavor
                        else "unknown"
                    ),
                    mass_definition=(
                        "effective electron-neutrino mass from tritium beta decay"
                        if is_electron_flavor
                        else "individual flavor mass not measured"
                    ),
                    source=_PDG_LEPTONS,
                    antiparticle=antiparticle,
                )
            )

    # The photon mass is experimentally bounded; the gluon masslessness is a
    # theoretical property of the Standard Model (the gluon has 8 color states).
    rows.extend(
        [
            _particle(
                "gamma",
                "photon",
                "gauge boson",
                0,
                mass_upper_limit_gev_c2=1.0e-27,
                mass_status="upper_limit",
                mass_definition="experimental upper limit",
                source=_PDG_BOSONS,
            ),
            {
                **_particle(
                    "g",
                    "gluon",
                    "gauge boson",
                    0,
                    mass_gev_c2=0.0,
                    mass_status="theoretical",
                    mass_definition="massless in the Standard Model; color octet",
                    source=_PDG_BOSONS,
                ),
                "state_count": 8,
            },
            _particle(
                "W+",
                "W+ boson",
                "gauge boson",
                1,
                80.3625,
                0.0077,
                source=_PDG_BOSONS,
            ),
            _particle(
                "W-",
                "W- boson",
                "gauge boson",
                -1,
                80.3625,
                0.0077,
                source=_PDG_BOSONS,
                antiparticle=True,
            ),
            _particle(
                "Z",
                "Z boson",
                "gauge boson",
                0,
                91.1879,
                0.0020,
                source=_PDG_BOSONS,
            ),
            _particle(
                "H",
                "Higgs boson",
                "scalar boson",
                0,
                125.13,
                0.11,
                source=_PDG_BOSONS,
            ),
        ]
    )

    dataframe = pd.DataFrame(rows)
    if additional_data is not None:
        if not isinstance(additional_data, Mapping):
            raise TypeError(
                "additional_data must map particle symbols to data mappings"
            )
        known_symbols = set(dataframe["symbol"])
        unknown_symbols = set(additional_data) - known_symbols
        if unknown_symbols:
            raise ValueError(f"Unknown particle symbols: {sorted(unknown_symbols)}")
        extra_rows = {}
        for symbol, values in additional_data.items():
            if not isinstance(values, Mapping):
                raise TypeError(f"Additional data for {symbol!r} must be a mapping")
            extra_rows[symbol] = dict(values)
        extra_dataframe = pd.DataFrame.from_dict(extra_rows, orient="index")
        extra_dataframe.index.name = "symbol"
        dataframe = dataframe.join(extra_dataframe, on="symbol", rsuffix="_additional")

    if source == False:
        dataframe = dataframe.drop("source", axis=1)

    return dataframe
