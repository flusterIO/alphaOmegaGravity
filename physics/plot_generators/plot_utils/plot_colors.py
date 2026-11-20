from dataclasses import dataclass, field
from typing import List
import random

"""Shared color palette for paper figures."""


@dataclass
class PlotColors:
    """CSS color strings used by the paper's plots.

    The first four colors are the existing Okabe-Ito-inspired paper series
    colors, in their established order. The remaining colors extend that
    palette for figures with more series.
    """

    colors: List[str] = field(
        default_factory=lambda: [
            "#0072B2",  # blue
            "#D55E00",  # vermilion
            "#009E73",  # bluish green
            "#E69F00",  # orange
            "#CC79A7",  # reddish purple
            "#56B4E9",  # sky blue
            "#F0E442",  # yellow
            "#222222",  # charcoal
        ]
    )

    def by_index(self, idx: int):
        return self.colors[idx % self.colors.__len__()]

    def random(self):
        return random.choice(self.colors)
