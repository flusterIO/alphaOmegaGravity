from dataclasses import dataclass
from typing import Union
import scipy.constants as spc
import numpy as np

[
    "name",
    "class",
    "mass_kg",
    "radius_m",
    "relative_velocity_m_s",
    "velocity_definition",
    "source",
    "compactness_GM_Rc2",
    "delta_s^-2",
    "omega_delta_cuberoot_s^-2_3",
    "source_url",
]


@dataclass
class MassiveBody:
    name: str
    _class: str
    # The mass in kilograms
    mass: float
    # The radius in meters
    radius: float

    relative_velocity: float
    velocity_definition: str
    rotation_period_s: Union[float, None]
    rotation_period_kind: str
    rotation_direction: str
    rotation_source: str
    rotation_source_url: str
    # The charge in coulumbs
    charge_value_C: Union[float, None]
    charge_value_type: Union[str, None]
    charge_abs_upper_bound_C: Union[float, None]
    charge_source: str
    charge_source_url: str
    charge_notes: str
    source: str
    source_url: str

    def compactness(self, G=spc.G, c=spc.c):
        return G * self.mass / (self.radius * c**2)

    def surface_g(self, G=spc.G):
        return G * self.mass / self.radius**2

    def delta(self, G=spc.G):
        return 2 * G * self.mass / self.radius**3

    def omega(self, G=spc.G):
        return (2 * G * self.mass / self.radius**3) ** (1 / 3)

    def derivedRelativeVelocity(self, G=spc.G):
        M = self.mass
        R = self.radius
        return np.sqrt(G * M * (G * M + 2 * R**3)) / (G * M + R**3)
