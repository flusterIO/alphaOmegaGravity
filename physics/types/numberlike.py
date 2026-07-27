from typing import TypeAlias, Union
import numpy.typing as npt


NumberLike: TypeAlias = Union[float, npt.NDArray]

OptionalFloat: TypeAlias = Union[float, None]
