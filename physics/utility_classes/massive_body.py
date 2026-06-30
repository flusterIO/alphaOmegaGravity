from dataclasses import dataclass
from physics.utility_classes.pos3 import Pos3


@dataclass
class MassiveBody:
    mass: float
    position: Pos3
