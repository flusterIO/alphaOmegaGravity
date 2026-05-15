from typing_extensions import Union
from physics.functions.statistics.statistic_utils import percentDiff
import astropy.units as u
from physics.symbolic.function_generators.theory.gravity.omega import omega
from physics.symbolic.function_generators.theory.quantum.alpha import alphaSolar
from physics.symbolic.function_generators.theory.velocity.vbar import (
    vbarCosmological,
    vbarKinematic,
)
import sympy as smp
from physics.symbolic.variables import subs, alpha, subsSolar
import scipy.constants as spc
from pathlib import Path
from physics.parameters import cmbRelativeVelocities


def describe_results(write_output_to: Union[str, None] = None):
    _vbarCosmological = vbarCosmological()
    _vbarKinematic = vbarKinematic()
    _omega = omega()
    _spinorApproximation = 1 / (2 - alpha)
    _alphaSolar = alphaSolar()
    omegaSubbed = _omega.subs(subs)
    spinorSubbed = _spinorApproximation.subs([[alpha, spc.alpha]])
    vbarKinematicSubbed = _vbarKinematic.subs(subs)
    alphaCorrection = alpha * (2 - alpha)
    s = f"""# Divergent Gravity 

This summary is generated at runtime, using CODATA and IAU values provided by the scipy & astropy packages.

### Defined Quantities

$$
\\omega = \\sqrt[3]{{\\frac{{d}}{{dR}} g}} = {smp.latex(_omega)} ~\\text{{m}} = {omegaSubbed:.8f} ~ \\text{{m}} = {smp.latex(alphaCorrection)} = {percentDiff(alphaCorrection.subs([[alpha, spc.alpha]]), omegaSubbed).evalf():.4f}\\%
$$

$$
\\Gamma = {smp.latex(_spinorApproximation)} = {spinorSubbed:.8f} ~\\frac{{1}}{{\\text{{m}}}}
$$

$$
\\Gamma \\omega = {(spinorSubbed * omegaSubbed).evalf()} = \\alpha \\pm {percentDiff((spinorSubbed * omegaSubbed).evalf(), spc.alpha):.4f}\\%
$$

### Velocities

The accuracy of this cosmological velocity is difficult to estimate due to the wildly differing results of observational research, but this falls almost perfectly on top of probability curves found through research correlating distant SNe to the microwave background.

$$
v_{{\\text{{cosmological}}}} = {smp.latex(_vbarCosmological)} = {_vbarCosmological.subs(subs):.2f} ~\\text{{m}}
$$

$$
v_{{\\text{{kinematic}}}} = {smp.latex(_vbarKinematic)} = {vbarKinematicSubbed:.2f}  ~\\text{{m}} = v_{{\\text{{observation}}}} \\pm {percentDiff(cmbRelativeVelocities[2], vbarKinematicSubbed * u.m / u.s):.2f}-{abs(percentDiff(cmbRelativeVelocities[0], vbarKinematicSubbed * u.m / u.s)):.2f}\\%
$$

$$
\\alpha_{{\\text{{Solar}}}} = {smp.latex(_alphaSolar)} = {_alphaSolar[0].subs(subsSolar).subs(subs):.8f}, {_alphaSolar[1].subs(subsSolar).subs(subs):.8f} = 1 \\pm \\alpha \\pm {percentDiff(spc.alpha + 1, _alphaSolar[0].subs(subsSolar).subs(subs)):.4f}\\%
$$

    """

    if write_output_to is None:
        return s
    else:
        p = Path(write_output_to)
        p.write_text(s)
        return s


describe_results(
    "/Users/bigsexy/Desktop/python/classes/physics/test_describe_output.md"
)
