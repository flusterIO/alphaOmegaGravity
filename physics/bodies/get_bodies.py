from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as spc
import pandas as pd
from dataclasses import dataclass, field
from physics.bodies.massive_body import MassiveBody
from physics.data.subsets import outliers


def get_bodies(drop_rows_by_name: List[str] = []):
    df = pd.read_csv(
        "/Users/bigsexy/Desktop/notes/content/resources/csv/gravity_geometry_dataset.csv"
    )
    bodies = []
    for i, body in df.iterrows():
        if not drop_rows_by_name.__contains__(body["name"]):
            massive_body = MassiveBody(
                name=body["name"],
                _class=body["class"],
                mass=body["mass_kg"],
                radius=body["radius_m"],
                relative_velocity=body["relative_velocity_m_s"],
                velocity_definition=body["velocity_definition"],
                source=body["source"],
                source_url=body["source_url"],
                rotation_period_s=(
                    body["rotation_period_s"] if "rotation_period_s" in body else None
                ),
                rotation_period_kind=(
                    body["rotation_period_kind"]
                    if "rotation_period_kind" in body
                    else None
                ),
                rotation_direction=(
                    body["rotation_direction"] if "rotation_direction" in body else None
                ),
                rotation_source=(
                    body["rotation_source"] if "rotation_source" in body else None
                ),
                rotation_source_url=(
                    body["rotation_source_url"]
                    if "rotation_source_url" in body
                    else None
                ),
                charge_value_C=(
                    body["charge_value_C"] if "charge_value_C" in body else None
                ),
                charge_value_type=(
                    body["charge_value_type"] if "charge_value_type" in body else None
                ),
                charge_abs_upper_bound_C=(
                    body["charge_abs_upper_bound_C"]
                    if "charge_abs_upper_bound_C" in body
                    else None
                ),
                charge_source=(
                    body["charge_source"] if "charge_source" in body else None
                ),
                charge_source_url=(
                    body["charge_source_url"] if "charge_source_url" in body else None
                ),
                charge_notes=body["charge_notes"] if "charge_notes" in body else None,
                # The charge in coulum=body[""],
            )
            bodies.append(massive_body)
    return bodies
