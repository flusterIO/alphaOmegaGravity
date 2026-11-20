from enum import Enum
from typing import Self
from IPython.core.getipython import get_ipython
from PIL import Image
import pyperclip as pc


class Images(Enum):
    CristoffelDiagram = "cristoffel-diagram"
    FiveDimensionalGeometry = "5d-geometry-diagram"
    EquitemporalPlaneGeometry = "equitemporal-plane-geometry"
    CovariantDerivativeVsCristoffel = "covariant-derivative-vs-cristoffel"

    def to_path(self):
        if self == Images.CristoffelDiagram:
            return "/Users/bigsexy/Desktop/swift/Fluster/apps/website/public/assets/diagrams/ai_generated/cristoffel_symbols.png"
        if self == Images.FiveDimensionalGeometry:
            return "/Users/bigsexy/Desktop/swift/Fluster/apps/website/public/assets/diagrams/ai_generated/5dGeometry.png"
        if self == Images.EquitemporalPlaneGeometry:
            return "/Users/bigsexy/Desktop/swift/Fluster/apps/website/public/assets/diagrams/ai_generated/equi_temporal_plane_geometry.png"
        if self == Images.CovariantDerivativeVsCristoffel:
            return "/Users/bigsexy/Desktop/swift/Fluster/apps/website/public/assets/diagrams/ai_generated/covariant_derivative_and_the_cristoffel_symbol.png"
        return "Invalid Path"

    def copy_path(self):
        pc.copy(self.to_path())


def renderImage(img: Images, open=False):
    p = img.to_path()
    ipy = get_ipython()
    if ipy is not None:
        ipy.run_line_magic("icat", p)
    if open:
        img = Image.open(p)  # pyright: ignore
        img.show()  # pyright: ignore
