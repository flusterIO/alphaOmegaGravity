import scipy.constants as spc
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Union, List

outliers = [
    "Milky Way",
    "Andromeda (M31)",
    "Typical neutron star",
    "Typical white dwarf",
]


def getDataFrame(drop_outliers: Union[bool, List[str]] = False):
    file_name = Path(__file__).parent / "gravity_geo_dataset.csv"
    df = pd.read_csv(file_name)
    if drop_outliers:
        drop_rows = outliers if drop_outliers == True else drop_outliers
        for r in drop_rows:
            df = df[(df["name"] != r)].copy()
    df["mass_density"] = df["mass_kg"] / (4 / 3 * spc.pi * df["radius_m"] ** 3)
    df["omega"] = (2 * spc.G * df["mass_kg"] / df["radius_m"] ** 3) ** (1 / 3)
    df["rotation_period_earth_day"] = df["rotation_period_s"].to_numpy() / np.array(
        df.iloc[[3]]["rotation_period_s"]
    )
    df["g_acceleration"] = spc.G * df["mass_kg"] / df["radius_m"] ** 2
    df["orbital_gamma"] = 1 / np.sqrt(
        1 - df["orbital_speed_mean_m_s"].to_numpy() ** 2 / spc.c**2
    )
    df["macro_alpha_equivalent"] = df["orbital_speed_mean_m_s"] / spc.c
    df["d_theta"] = 2 * spc.pi / df["rotation_period_s"]
    return df
