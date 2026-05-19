from typing_extensions import Union
import astropy.constants.iau2015 as ac
from physics.functions.statistics.statistic_utils import percentDiff
import astropy.units as u
from physics.parameters import earthRadiusPolar, earthRadiusMean
from physics.symbolic.function_generators.theory.em.current_equivalent import currentEquivalent
from physics.symbolic.function_generators.theory.em.electric_field_divergence import (
    electricFieldDivergence,
)
from physics.symbolic.function_generators.theory.gravity.delta import delta
from physics.symbolic.function_generators.theory.gravity.mass_density import (
    massDensityEquivalence,
)
from physics.symbolic.function_generators.theory.gravity.omega import omega
from physics.symbolic.function_generators.theory.quantum.alpha import alphaSolar
from physics.symbolic.function_generators.theory.velocity.vbar import (
    vbarAlphaEquivalence,
    vbarCosmological,
    vbarKinematic,
)
import sympy as smp
from physics.symbolic.variables import subs, alpha, subsSolar
import scipy.constants as spc
from pathlib import Path
from physics.parameters import cmbRelativeVelocities


def describe_results(write_output_to: Union[str, None] = None):
    """
    Prints a markdown file, templated at runtime using CODATA and IAU constants
    provided by the scipy and astropy package.
    """
    _vbarCosmological = vbarCosmological()
    _vbarKinematic = vbarKinematic()
    _omega = omega()
    _spinorApproximation = 1 / (2 - alpha)
    _alphaSolar = alphaSolar()
    omegaSubbed = _omega.subs(subs)
    spinorSubbed = _spinorApproximation.subs([[alpha, spc.alpha]])
    vbarKinematicSubbed = _vbarKinematic.subs(subs)
    alphaCorrection = alpha * (2 - alpha)
    vbarAlpha = vbarAlphaEquivalence()
    electricDivergence = electricFieldDivergence()
    electricDivergenceSubbed = electricDivergence.subs(subs)
    rotV = 2 * spc.pi * earthRadiusPolar / 86400
    md = massDensityEquivalence()
    rhoSubbed = md[0].subs(subs).subs([[alpha, spc.alpha]]).evalf()
    _delta = delta()
    vbarAlphaDiff = percentDiff(vbarAlpha[0].subs(subs), (spc.alpha * (2 - spc.alpha))**3)
    vbarKinematicDiff = percentDiff(cmbRelativeVelocities[2], vbarKinematicSubbed * u.m / u.s)
    vbarKinematicDiff1 = percentDiff(cmbRelativeVelocities[0], vbarKinematicSubbed * u.m / u.s)
    rhoDiff = percentDiff(ac.M_earth.value.tolist() / (4/3 * spc.pi * earthRadiusMean**3), rhoSubbed)
    jDiff = abs(percentDiff(_delta.subs(subs) * _vbarCosmological.subs(subs), spc.golden_ratio))
    currentEquivDiff = abs(percentDiff(currentEquivalent().subs(subs), spc.golden_ratio))
    s = f"""# Divergent Gravity 

This is a short summary generated at runtime, using CODATA and IAU values for all equations provided by the scipy & astropy packages.

There is a short derivation included, but for a longer derivation, please see the documentation at [flusterapp.com](https://flusterapp.com/docs). There is a more complete derivation included as an example for the application.


| Derivation                                                | Function                         | Error                              |
| --------------------------------------------------------- | -------------------------------- | ---------------------------------- |
| cosmological velocity from local gravitational parameters | `vbarCosmological`               | Difficult to assess                |
| kinematic velocity from local gravitational parameters    | `vbarKinematic`                  | ${vbarKinematicDiff:.3f}\\% - {vbarKinematicDiff1:.3f}\\%$                  |
| Earth's rotational velocity via $\\alpha$                  | `rotationVelocityFromDivergence` | $0.006778\\%$                       |
| $\\alpha$, $2-\\alpha$, $1 \\pm \\alpha$                      | `alpha`                          | $0.025\\%$ <br/> $0.0034\\%$ (solar) |
| Earth's cosmological velocity as a function of $\\alpha$   | `vbarAlphaEquivalence`           | ${vbarAlphaDiff:.5f}\\%$                           |
| $\\Phi$ mass density approximation   | `massDensityEquivalence`           | ${rhoDiff:.5f}\\%$                           |
| $\\mathbf{{J}} = \\Phi$ (A radius exists on Earth's surface where both are equivalent)   | `currentEquivalent`           | ${currentEquivDiff:.5f}\\%$                           |


### The Lighthouse & The Clocktower

> A quick thought experiment to justify a case for classical synchronicity

Consider two observers, $A$ and $B$. Let $A$ and $B$ agree that $B$ should travel between two arbitrary points in space (assuming a uniform gravitational field), $p_1$ and $p_2$ at an _agreed_ upon velocity $v$. 

At $p_2$ place a time keeping device that remains visible from $p_1$. Allow that $A$ remains at rest relative to $p_1$, and very near $p_1$ while $B$ goes in motion (at an inertial pace) between $p_1$ and $p_2$. 

At what time does $A$ observe $B$ reach $p_2$? If the value is not $d/v$, but $\\gamma d/v$, then $v \\neq v_\\prime$. The only possible solution that allows for Einstein's additional factor of $\\gamma$ while allowing $A$ to remain in an unchanged state is that $d \\to d\\gamma$ in the coordinate system of the observer in relative motion, while $v \\to v/\\gamma$ in the coordinate system of the observer at relative rest to accomodate the change in spatial density.


---

## Kinematic & Cosmological Velocities

From this simple thought experiment, it's trivial to derive our cosmological velocity... simple solve $\\gamma$ for $v$ and plugin in the average of the integral of $\\delta$, the derivative of gravitational acceleration with respect to radius:

$$
\\delta = \\frac{{d}}{{dR}} g = 2G \\frac{{M}}{{R^3}}
$$

This gives:

$$
v_{{\\text{{cosmological}}}} = {smp.latex(_vbarCosmological)} = {_vbarCosmological.subs(subs):.2f} ~\\text{{m}}
$$

The accuracy of this cosmological velocity is difficult to estimate due to the wildly differing results of observational research, but this falls almost perfectly on top of probability curves found through research correlating distant SNe to the microwave background. I'll cite work properly when I have access to the internet reliably again and can dedicate the time to it, but see the work of A. Singal, K. Land, or C. Gordon.

From here, by using 'time-equivalent' quantities... ratios of distance/distance of a magnitude equivalent to $t$ we can derive the following for our kinematic velocity:

$$
v_{{\\text{{kinematic}}}} = {smp.latex(_vbarKinematic)} = {vbarKinematicSubbed:.2f}  ~\\text{{m}} = v_{{\\text{{observation}}}} \\pm {percentDiff(cmbRelativeVelocities[2], vbarKinematicSubbed * u.m / u.s):.2f}-{abs(percentDiff(cmbRelativeVelocities[0], vbarKinematicSubbed * u.m / u.s)):.2f}\\%
$$

As I'm without internet I'm unable to look up the paper who's velocity estimates I'm using, but it is one of the most cited papers producing related data of the past several years.

---

## Electromagnetism

If this model is correct, and Earth dilates _with_ the rest of the Universe, along with all other bodies, than we should consider the divergence of $\\delta$... or rather, $\\delta$ **as** divergence.

Taking the inverse of $\\nabla \\cdot \\omega = \\delta$, assuming a diagonal 3-dimensional matrix gives:

$$
\\omega = \\sqrt[3]{{\\frac{{d}}{{dR}} g}} = {smp.latex(_omega)} ~\\text{{m}} = {omegaSubbed:.8f} ~ \\text{{m}} = {smp.latex(alphaCorrection)} = {percentDiff(alphaCorrection.subs([[alpha, spc.alpha]]), omegaSubbed).evalf():.4f}\\%
$$

This gives the equation that I've spent the past 5 years looking for:

$$\\boxed{{
\\alpha = \\frac{{1}}{{2 - \\alpha}} \\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}}
}}
$$

Yes, that's $\\alpha$, the fine-structure constant, but don't worry... we can do better. Now that we know it's here and how to treat the math thanks to the work already done on electromagnetism, it pops up **everywhere**.

And no, this is **not** dimensionally incorrect. A fundamental part of this model is that the temporal axis _is_ cosmic inflation. By combining these two processes into a single process, we must remove a time derivative and replace the fundamental, driving rate of change with motion... not time. This then takes the dimensionality of:

$$
\\overbrace{{\\alpha}}^{{dt}} = \\underbrace{{\\frac{{1}}{{2 - \\alpha}}}}_{{1/dx}} \\overbrace{{\\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}}}}^{{ds}}
$$

Consider the function for our Sun, where our 'Spinor' approximation limits at:

$$
\\lim_{{\\Gamma \\to (1/2)^{{1/3}}}} \\frac{{1}}{{2 - \\alpha}} \\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}} \\approx 1 \\pm \\alpha
$$

Or, from Python and Astropy:

$$
\\alpha_{{\\text{{Solar}}}} = {smp.latex(_alphaSolar)} = {_alphaSolar[0].subs(subsSolar).subs(subs):.8f}, {_alphaSolar[1].subs(subsSolar).subs(subs):.8f} = 1 \\pm \\alpha \\pm {percentDiff(spc.alpha + 1, _alphaSolar[0].subs(subsSolar).subs(subs)):.4f}\\%
$$

> Ignore the goofy exponents if they appear astronomical. That's coming directly from sympy, but they evaluate to 1/3.

Above, 'dt' is used simply to explain the phenomena in terms we are familiar with. There is absolutely no need for a temporal axis, and it plays absolutely no role in the underlying physics. It is simply a dimensionless ratio of change in the form of $ds/dx$, giving rise to the **sensation** of time due to our comotion and the dilation of $c$ along this axis. There is no drag, no resistance or boost... no physical interaction with the 'time medium' of any kind. There is no such thing.

Riddle me this... how can $c$ be constant in-between reference frames, when it's a time dependent quantity and time dilates? This logic is self defeating. _If_ $c$ is of an equivalent magnitude ***in*** each reference frame, than it ***must*** dilate between reference frames to accomodate the shift in unit magnitude. It cannot possibly be both, and this is why $c$ _must_ dilate between reference frames. It's inherent in the components of $c$, $\\mu_0$ and $\\epsilon_0$... the permitivity and permability of _free space_, but with units proportional to $1/\\text{{length}}$ Of course they're going to respond to a shift in spatial density.


Let's lay this out more explicitly, since many pseudo-academics remain determined to pray to their volcano.

Let's consider the divergence of the magnetic field. If this model is correct, than the divergence of the magnetic field $\\neq 0$, but rather $\\delta$.

Assuming the identity:

$$
\\nabla \\cdot E = c \\left(\\nabla \\cdot B \\right)
$$

This gives:

$$
{smp.latex(electricDivergence)} = {electricDivergence.subs(subs)}
$$

Using our 'spinor' approximation $\\Gamma$, we find:


$$
\\Gamma = {smp.latex(_spinorApproximation)} = {spinorSubbed:.8f} ~\\frac{{1}}{{\\text{{m}}}}
$$
Where:

$$
\\Gamma \\omega = {(spinorSubbed * omegaSubbed).evalf()} = \\alpha \\pm {percentDiff((spinorSubbed * omegaSubbed).evalf(), spc.alpha):.4f}\\%
$$

$$
{smp.latex(electricDivergence)} = {spinorSubbed * electricDivergenceSubbed} ~\\text{{m}} ~\\text{{s}}^{{-1}}
$$

This is our rotational velocity to within:

$$
\\frac{{2 \\pi R_\\oplus}}{{86400}} \\pm {percentDiff(rotV, spinorSubbed * electricDivergenceSubbed):.4f}\\%
$$

And a direct, mechanistic link to why our time dilates at the poles due to gravity at an almost identical rate to that of the equator due to the additional rotational velocity. Earth is *not* special, but these units are *not* constant.

This is off by ${rotV - (spinorSubbed * electricDivergenceSubbed):.5f}$ meters per second... well within measurement error.

> Note here that our polar radius was used to avoid the influence of centripital force.

Then, consider the electromagnetic displacement current $\\mathbf{{J}}$ where

$$
\\mathbf{{J}} = \\rho_e v_0
$$

Since, in electromagnetism:

$$
\\nabla \\cdot E = \\frac{{\\rho_e}}{{\\epsilon_0}}
$$

Then

$$
\\mathbf{{J}} \\propto \\delta \\cdot v_0
$$

Or

$$
\\mathbf{{J}} \\propto {_delta.subs(subs) * _vbarCosmological.subs(subs)}
$$

or 

$$
\\mathbf{{J}} \\propto \\Phi \\pm {jDiff:.4f}\\%
$$

---

Finally, consider the relationship:

$$
v_0^2 = {smp.latex(vbarAlpha[0])} = {smp.latex(vbarAlpha[1])} \\pm {vbarAlphaDiff:.5f}
$$

Then let's consider that

$$
\\rho = \\frac{{M}}{{V}} = \\frac{{M}}{{\\frac{{4}}{{3}} \\pi R^3}} = \\frac{{3M}}{{4\\pi R^3}}
$$

And the fact that we defined 

$$
\\delta = 2 G \\frac{{M}}{{R^3}}
$$

Then

$$
\\frac{{M}}{{R^3}} = \\frac{{\\delta}}{{2G}}
$$

and

$$
\\rho = \\frac{{3}}{{4\\pi}} \\left(\\frac{{\\delta}}{{2G}}\\right) = \\frac{{3 \\delta}}{{8 \\pi G}}
$$

And given the relationship that

$$
\\delta \\approx \\frac{{\\Phi}}{{v_0}}
$$

then finally


$$
\\rho = \\frac{{3 \\left( \\frac{{\\delta R^{{3}}}}{{2G}} \\right)}}{{4 \\pi R^{{3}}}} = \\frac{{3 \\delta}}{{8 \\pi G}}
$$

Since we found that the dimensionless ratio $ds/dx = \\delta v_{{0}} \\approx \\Phi$, then

$$
\\rho = \\frac{{\\Phi}}{{\\left( \\frac{{8\\pi}}{{3}}G \\cdot v_{{0}} \\right)}}
$$

Using the constraint we found above, that $\\beta = \\left( \\alpha \\left( 2 - \\alpha \\right) \\right)^{{2/3}}$

$$
\\rho = \\frac{{\\delta}}{{G}}\\cdot \\frac{{\\Phi/\\left( c \\delta \\right)}}{{\\alpha \\left( 2 - \\alpha \\right)}}
$$

Or

$$
\\rho = \\frac{{\\Phi}}{{Gc \\cdot \\alpha \\left( 2 - \\alpha \\right)}}
$$

This is our traditionally calculated mean density, $\\rho_\\prime$

$$
\\rho_\\prime \\pm {rhoDiff:.5f}\\%
$$

Not a _great_ match, but pretty good considering almost every component in there is an approximation, slid by some $\\Delta$. My point is not to say that our mass density can be found by this function, but rather that the geometry of our Universe can be better understood by studying these approximations. _One_ of these alignments may be a 'dimensional coincidence', as one email put it... but not **all** of them.

And no, the dimensionality is not off. There is no such thing as time, and we need to remove it from our text books. Our 4th displacement coordinate is best represented by a position along the density axis that arises from this geometry, and the _driving_ rate of change is and always has been _motion..._ not time. Time is nothing but a useful mathematical tool.


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
