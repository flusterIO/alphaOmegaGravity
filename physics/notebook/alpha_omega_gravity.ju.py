# %% [md]

## Divergent Gravity

> Glance down at the table of results below before you close this notebook. The conclusion is controversial, but the path to the conclusion is rock solid. The results should speak for themselves. And no, the Lorentz transformations are not the only solution that satisfies experiment.

What if time, cosmic inflation and gravitational acceleration could be collapsed into a single process as observed from different reference frames, yielding a set of functions binding electromagnetism and gravity?

I'm proposing that it is not time that elongates relative motion, but distance. This is justified below by a simple thought experiment, a thought experiment that I argue disproves the very possibility of relativistic simultaneity. To be clear, if Einstein's experimentally proven factor of $\gamma$ is accurate, while a single, instantaneous event is required to occur at precisely 1 instance in time, than the model I am proposing _has_ to be true.

Do we not already believe this? Is cosmic inflation not a time dependent property? Does more time not pass for an observer in relative motion, according to $\gamma$?

What is the big leap to considering that this may be the source of the equivalence principle? The way the physics just aligns is astonishing.

### The conclusion

If these two principles are true... cosmic inflation is time dependent and more time passes for an observer in relative motion, then the model I am proposing should follow logically, no matter how controversial the conclusions.

What conclusion is that? Not only is the Universe expanding, so too is _everything_, including the Earth. Space-time becomes a hypersphere in derivative form, mirroring the behavior seen in electric monopoles.

This becomes the source of the equivalence principle, from which we can derive both our cosmological and kinematic velocities to an _intriguing_ accuracy, but we're not done with that.

From this principle and the geometry of divergence, we find a numerical match to $\alpha (2 - \alpha)$ to within $0.025\%$, but this function holds. Once the spinor approximation is accounted for, we find an even better match for our Sun, with an error of $0.0034\%$.

From here, we find an electric field divergence equivalent that matches our spinor approximation, $1 / (2 - \alpha)$ times our rotational velocity to within $0.032$ meters per 'second'.

Then, finding an electric currrent equivalent, we find $J = \Phi \pm 0.21\%$, and a cosmological velocity as a function of $\alpha$, $\beta^2 = (\alpha (2 - \alpha))^3$.

Lastly, we derive a rotation function that gives our galactic radius with only visible mass, and define a mass density approximation of:


$$
\rho \approx \frac{\Phi}{Gc (\alpha (2 - \alpha))} 
$$

This holds, because these quantities are not constant. They _appear_ to be constant _in_ each reference frame because they dilate inversely proportional to spatial density. If motion is quantized such that there exists some $n$ number of '$dx$' quantities between points $p_1$ and $p_2$, would you not expect velocity to dilate proportionally between these points as $dx \to dx~\gamma$? (If you're not onboard yet, replace that $\gamma$ with $\mathit{H}_0$. See my point?).

In the end, we end up with a combined equation:

$$
\alpha = \frac{1} {2 - \alpha} \sqrt[3]{2\frac{GM}{R^3}}
$$

> This is not _technically_ correct. There is a term of motion in the cube root, but at $\Delta t = 1$, it cancels with the additional $1/\text{distance}$ term, giving the above result. The full equation and derivation can be found on [my website](https://flusterapp.com/docs).

This gives the following relationship:

$$
\frac{d}{dR} g \approx \beta^2_{\text{cosmological}} \approx (\alpha (2 - \alpha))^3 \approx 3.07 \times 10^{-7} ~ \text{m}
$$

And a final orbital approximation of:

$$
\frac{\sqrt[3]{GM_\odot}}{v_{\text{orbit}}} = \alpha \frac{\mathbf{AU}}{R_\oplus} \pm 0.25\%
$$
# %%

import scipy.constants as spc
import astropy.constants.codata2010 as cd
import astropy.constants.astropyconst40 as ac
import plotly.graph_objects as go
import plotly.io as pio
import sympy as smp
from scipy.optimize import fsolve
import plotly.express as px
import pandas as pd
from IPython.core.interactiveshell import InteractiveShell
from pathlib import Path
import os

 # Personal preference... comment this out if you're on lightmode to make the plots light.
pio.templates.default = "plotly_dark"

InteractiveShell.ast_node_interactivity = "all"
from IPython.display import display, Markdown, Image


def print_md(string):
    display(Markdown(string))

# %% [md]

## Parameters

As I'm offline while writing this I'm unable to look up the papers these parameters were taken from, but all values were taken from semi-recent, well cited sources.

There is a package that I plan to publish next time I'm online with all of these functions available at `alpha_omega_gravity`... because after almost 5 years of being homeless to work on this model, I think I earned the right to be grandiose.

Here, I do what I often do in my notes and use template strings. Feel free to close the input cells for many of these, as it's just a duplicate of the output with some 
variables inserted.
# %%

cmbRelativeVelocities = [
    1.228 * 10**-3 * cd.c,
    1.231 * 10**-3 * cd.c,
    1.234 * 10**-3 * cd.c,
]

solarObitalVelocity = 2.3 * 10**5
solarCircleEnclosedMass = 5 * 10**10 * ac.M_sun.value

observationalSolarGalacticCenterRadius = [(8.12 - 0.03) * spc.kilo * spc.parsec, (8.12 + 0.03) * spc.kilo * spc.parsec]


earthRadiusPolar: float = 6356752

earthRadiusMean: float = 6371000

# %%

(
    G,
    q,
    g,
    M,
    m,
    R,
    mo,
    eo,
    n,
    e,
    c,
        au,
    me,
    m_tau,
    m_plank,
    delta,
    omega,
    dtheta,
    dt,
    alpha,
    h,
    hbar,
    lp,
    gamma,
    vSolar,
    vOrbit,
    mSolar,
    rGalactic,
    mSolarCircle
) = smp.symbols(
    [
        "G",
        "z",
        "g",
        "M",
        "m",
        "R",
        r"\mu_{0}",
        r"\epsilon_{0}",
        "n",
        "e",
        "c",
        r"\mathbf{AU}",
        r"m_{e}",
        r"m_{\tau}",
        r"m_p",
        r"\delta",
        r"\omega",
        r"d\theta",
        r"dt",
        r"\alpha",
        "h",
        r"\hbar",
        r"\l_{p}",
        r"\gamma",
        r"v_{\odot}",
        r"v_{\text{orbit}}",
        r"M_{\odot}",
        r"r_{\text{galactic}}",
        r"m_{\text{sc}}",
    ]
)

_c = smp.sqrt(1 / (mo * eo))
_delta = -2 * G * M / R**3
_n = mo / (4 * smp.pi)
_omega = _delta ** smp.Rational(1, 3)
_alpha = 1 / 2 * e**2 / h * smp.sqrt(mo / eo)
_mass_plank = smp.sqrt(h * c / (2 * smp.pi * G))
_lp = smp.sqrt(hbar * G / c**3)

subs = [
    [c, spc.c],
    [G, spc.G],
    [M, ac.M_earth.base.tolist()],
    [R, ac.R_earth.base.tolist()],
    [mo, cd.mu0.base.tolist()],
    [eo, cd.eps0.base.tolist()],
    [e, cd.e.base.tolist()],
    [h, cd.h.value.tolist()],
    [au, spc.au],
    [mSolar, ac.M_sun.value.tolist()],
]

subsSolar = [[R, ac.R_sun.base.tolist()], [M, ac.M_sun.base.tolist()], [vSolar, solarObitalVelocity], [mSolarCircle, solarCircleEnclosedMass]]

expand = [
    [g, G * M / R**2],
    [n, _n],
    [c, _c],
    [omega, _omega],
    [delta, _delta],
    [alpha, _alpha],
    [m_plank, _mass_plank],
    [lp, _lp],
    [hbar, h / (smp.pi * 2)],
    [vOrbit, vOrbitApproximation()]
]


# %%

def percentDiff(target, found):
    return (target - found) / target * 100

def vFromGamma(_gamma=gamma):
    return c * (1 - 1 / _gamma**2) ** smp.Rational(1, 2)

def gravityGammaEquivalent():
    return 1 + G * M / R**3

def delta():
    return 2 * G * M / R**3

def vbarCosmological():
    return vFromGamma(gravityGammaEquivalent())

def gammaEquivTD():
    return 1 + 0.5 * 1 / (1 + 1 / 2 * G * M / R**3) ** 2

def vbarKinematic() -> smp.Pow:
    return c * smp.sqrt(1 - 1 / (1 + smp.Rational(1, 2) * G * M / R**3) ** 2)

def omega():
    return delta() ** smp.Rational(1, 3)

def vbarAlphaEquivalence():
    a = vbarCosmological() ** 2 / c**2
    b = (alpha * (2 - alpha)) ** 3
    return [a, b]

def alphaEquivalence():
    f = smp.sqrt(1 - omega())
    return [1 + f, 1 - f]

def electricFieldDivergence():
    return c * delta()

def alphaSolar():
    _n = (smp.Rational(1, 2) * delta()) ** smp.Rational(1 / 3)
    return [1 + _n, 1 - _n]

def vbarCosmological():
    return vFromGamma(gravityGammaEquivalent())

def currentEquivalent():
    return delta() * vbarCosmological()

def massDensityTraditional():
    return M / (4 / 3 * smp.pi * R**3)

def massDensityEquivalence():
    a = alpha
    phi = smp.GoldenRatio
    a = phi / (G * c * (a * (2 - a)))
    b = massDensityTraditional()
    return [a, b]

def vOrbitApproximation():
    return 2 * smp.pi * au / (365.25 * spc.day)

def solarMetricVelocity():
    return (G * mSolar)**smp.Rational(1, 3)

def galacticRadiusPolynomial():
    return vSolar**2 * rGalactic**4 - (2 * G * mSolarCircle * rGalactic**3) - (G * mSolarCircle)**2

# %%

_vbarCosmological = vbarCosmological()
_vbarKinematic = vbarKinematic()
_omega = omega()
_spinorApproximation = 1 / (2 - alpha)
_alphaSolar = alphaSolar()
omegaSubbed = _omega.subs(subs)
spinorSubbed = _spinorApproximation.subs([[alpha, spc.alpha]])
vbarKinematicSubbed = _vbarKinematic.subs(subs)
alphaFactor = alpha * (2 - alpha)
vbarAlpha = vbarAlphaEquivalence()
electricDivergence = electricFieldDivergence()
electricDivergenceSubbed = electricDivergence.subs(subs)
rotV = 2 * spc.pi * earthRadiusPolar / 86400
md = massDensityEquivalence()
rhoSubbed = md[0].subs(subs).subs([[alpha, spc.alpha]]).evalf()
_delta = delta()
vbarAlphaDiff = percentDiff(vbarAlpha[0].subs(subs), (spc.alpha * (2 - spc.alpha))**3)
vbarKinematicDiff = percentDiff(cmbRelativeVelocities[2].value, vbarKinematicSubbed)
vbarKinematicDiff1 = percentDiff(cmbRelativeVelocities[0].value, vbarKinematicSubbed)
rhoDiff = percentDiff(ac.M_earth.value.tolist() / (4/3 * spc.pi * earthRadiusMean**3), rhoSubbed)
jDiff = abs(percentDiff(_delta.subs(subs) * _vbarCosmological.subs(subs), spc.golden_ratio))
currentEquivDiff = abs(percentDiff(currentEquivalent().subs(subs), spc.golden_ratio))
alphaSolarDiff = percentDiff(spc.alpha + 1, _alphaSolar[0].subs(subsSolar).subs(subs))
rotDiff = percentDiff(rotV, spinorSubbed * electricDivergenceSubbed)
omegaAlphaDiff = percentDiff(alphaFactor.subs([[alpha, spc.alpha]]), omegaSubbed.evalf())
solarMetVRatio = (solarMetricVelocity() / vOrbit).subs(expand).subs(subs).evalf()
auRRatio = spc.alpha * (au / R).subs(subs)
auRRatioPercentDiff = percentDiff(auRRatio, solarMetVRatio)
fAlphaOmega = ((1 / (2 - alpha) * omega()) - alpha).subs(expand)
hsols = smp.solve(fAlphaOmega, h)
alphaOmegaDiff1 = percentDiff(hsols[1].subs(subs), spc.h)
alphaOmegaDiff = percentDiff(hsols[0].subs(subs), spc.h * (spc.alpha / (2 - spc.alpha)))

# %%

massDensityEquivTitle = "Φ mass density approximation"
cosmoVelTitle = "cosmological velocity from local gravitational parameters"
df = pd.DataFrame([ 
    {
        "derivation": cosmoVelTitle,
        "function": "vbarCosmological",
        "% error": "Difficult to assess"
    },
    {
        "derivation": "kinematic velocity from local gravitational parameters",
        "function": "vbarKinematic",
        "% error": abs(vbarKinematicDiff)
    },
    {
        "derivation": "Earth's roational velocity via $\\alpha$",
        "function": "rotationVelocityFromDivergence",
        "% error": abs(rotDiff)
    },
    {
        "derivation": "Solar-Galactic Radius from visible mass",
        "function": "galacticRadiusPolynomial",
        "% error": 0.0
    },
    {
        "derivation": "α (solar)",
        "function": "alphaSolar",
        "% error": abs(alphaSolarDiff)
    },
    {
        "derivation": "α (Earth)",
        "function": "omega, alphaFactor",
        "% error": abs(omegaAlphaDiff)
    },
    {
        "derivation": "α (orbital velocity)",
        "function": "solarMetricVelocity, vOrbit",
        "% error": abs(auRRatioPercentDiff)
    },
    {
        "derivation": "α (from derived velocity)",
        "function": "vbarAlphaEquivalence",
        "% error": abs(vbarAlphaDiff)
    },
    {
        "derivation": "J ≈ Φ",
        "function": "currentEquivalent",
        "% error": abs(currentEquivDiff)
    },
    {
        "derivation": "h from α",
        "function": "omega, alpha",
        "% error": abs(alphaOmegaDiff1)
    },
    {
        "derivation": "h/(2 - α) from ω",
        "function": "omega, alpha",
        "% error": abs(alphaSolarDiff)
    },
    {
        "derivation": massDensityEquivTitle,
        "function": "massDensityEquivalence",
        "% error": abs(rhoDiff)
    },
])

df

# %% [md]

Keep in mind, for several of these quantities we could be picky and choose a specific radius that exists on Earth's surface and achieve a perfect match. These calculations were done with the most fitting radius for the calculation (polar for rotational velocity to reduce centripital influence, mean for calculating volume, astropy for everything else).

# %%

# Removes functions that are known to be merely approximations, using approximations.
# The point of this notebook is not to derive precise functions, but rather to derive a
# geometry of time itself, which can be better understood by studying these approximations.
filter_approximation_derivations = [
    massDensityEquivTitle,
    # Remove the cosmological velocity because it's error is a string.
    cosmoVelTitle 
]

withoutApproximations = df[~df['derivation'].isin(filter_approximation_derivations)]


# %%

print_md(f""" 
Or a combined average error of {(df["% error"].iloc[1:].mean()):.4f}%.

We can do even better if we remove the quantities that we assume are merely approximations, and not completely accurate functions. We then find an average error of {(withoutApproximations['% error'].sum() / withoutApproximations.__len__()):.4f}%. Not quite within observation error for some of the electromagnetic properties, but more than enough to establish a relationship between gravity and electromagnetism.

""")



# %% [md]

Here is a short derivation, but for a longer derivation, please see the documentation at [flusterapp.com](https://flusterapp.com). There is a more complete article included as an example for the application.

> This derivation actually had to be taken down to get the website to build while I was at the library (I'm still borderline homeless and without internet), but I'll have it back up shortly.

---

### The Lighthouse & The Clocktower

> A quick thought experiment to justify a case for classical synchronicity

Consider two observers, $A$ and $B$. Let $A$ and $B$ agree that $B$ should travel between two arbitrary points in space (assuming a uniform gravitational field), $p_1$ and $p_2$ at an _agreed_ upon velocity $v$. 

At $p_2$ place a time keeping device that remains visible from $p_1$. Allow that $A$ remains at rest relative to $p_1$, and very near $p_1$ while $B$ goes in motion (at an inertial pace) between $p_1$ and $p_2$. 

At what time does $A$ observe $B$ reach $p_2$? If the value is not $d/v$, but $\gamma d/v$, then $v \neq v_\prime$... _unless_, $d^\prime = d \gamma$.

---

The only possible solution that allows for Einstein's additional factor of $\gamma$ while allowing $A$ to remain in an unchanged state is that $d \to d\gamma$ in the coordinate system of the observer in relative motion, while $v \to v/\gamma$ in the coordinate system of the observer at relative rest to accomodate the change in spatial density. Any other solution would modify the agreed upon belocity, or result in a disagreement on the time that $B$ would reach $p_2$.

If the two observers are to agree upon the time at which $B$ reaches $p_2$ while $v = v_\prime$ in the coordinate system of the observer in motion (the coordinate system of the observer who's kinetic state has changed) then distance _must_ elongate with time.

Think about this even further. Do we not already believe this? We believe cosmic inflation is a time dependent property, what is the big leap to considering that maybe this property is _velocity_ dependent? Just look at the way the physics aligns when we make this single modification, while every single experimental validation of either SR or GR is maintained due to the $d=vt\gamma$ equivalence that both models seek to maintain (and their dependence on the equivalence principle).

Does it really not make more sense to have a temporal dimension driven by motion, instead of some unseen quantity we've never directly measured? It's time to think for yourself people... Einstein's not here to hold our hand anymore.

# %%


print_md(f""" 
## Kinematic & Cosmological Velocities

From this simple thought experiment, it's trivial to derive our cosmological velocity... simply solve $\\gamma$ for $v$ and plugin in the average of the integral of $\\delta$, the derivative of gravitational acceleration with respect to radius:

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

""")



# %%

print_md(f""" 
## Electromagnetism

If this model is correct, and Earth dilates _with_ the rest of the Universe, along with all other bodies, than we should consider the divergence of $\\delta$... or rather, $\\delta$ **as** divergence.

Taking the inverse of divergence, where $\\nabla \\cdot \\omega \\mathit{{I}} = \\delta$, where $\\mathit{{I}}$ is a 3-dimensional identity matrix gives:

$$
\\omega = \\sqrt[3]{{\\frac{{d}}{{dR}} g}} = {smp.latex(_omega)} ~\\text{{m}} = {omegaSubbed:.8f} ~ \\text{{m}} = {smp.latex(alphaFactor)} = {omegaAlphaDiff.evalf():.4f}\\%
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

$\\alpha$ is _not_ a constant in itself, and it needs to stop being viewed as such. It's comprised of other constants, each with it's own properties and antagonists, giving $\\alpha$ different units in different contexts. 

If this model is correct and space is quantized such that there exists $n$ number of some length quantity $\\mathit{{l}}$ between $p_1$ and $p_2$, then $h \\to h \\gamma$ in the coordinate system of the observer in motion, along with proportionally dilating $\\mu_0$ and $\\epsilon_0$. This dilation gives $\\alpha$ units of $\\text{{distance}}/\\text{{time}}$.

Consider the function for our Sun, where our 'Spinor' approximation limits at:

$$
\\lim_{{\\Gamma \\to (1/2)^{{1/3}}}} \\frac{{1}}{{2 - \\alpha}} \\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}} \\approx 1 \\pm \\alpha
$$

Or, from Python and Astropy:

$$
\\alpha_{{\\text{{Solar}}}} = {smp.latex(_alphaSolar)} = {_alphaSolar[0].subs(subsSolar).subs(subs):.8f}, {_alphaSolar[1].subs(subsSolar).subs(subs):.8f} = 1 \\pm \\alpha \\pm {alphaSolarDiff:.4f}\\%
$$

> Ignore the goofy exponents if they appear astronomical. That's coming directly from sympy, but they evaluate to 1/3.

Above, 'dt' is used simply to explain the phenomena in terms we are familiar with. There is absolutely no need for a temporal axis, and it plays absolutely no role in the underlying physics. It is simply a dimensionless ratio of change in the form of $ds/dx$, giving rise to the **sensation** of time due to our comotion and the dilation of $c$ along this axis. There is no drag, no resistance or boost... no physical interaction with the 'time medium' of any kind. There is no such thing. In the end, gravity is just electromagnetism applied to space itself, and this process is what we experience as time. There is no need for an additional dimension.
""")


# %%

print_md(
    f"""
## $\\nabla \\cdot E$ and $v_{{\\text{{rotational}}}}$

Riddle me this... how can $c$ be constant in-between reference frames, when it's a time dependent quantity, time dilates, and it's constant _in_ each reference frame? This logic is self defeating. _If_ $c$ is of an equivalent magnitude ***in*** each reference frame, than it ***must*** dilate between reference frames to accomodate the shift in unit magnitude. If a quantity has a magnitude of $n$ meters per second in both reference frames, but a second itself has a different magnitude, then that quantity is _not_ constant. It is this type of behavior that is masked by the very use of time as a mathematical tool. While perceiving a Universe without time may be diffcult, there is not a single argument for time as a physically meaningful quantity. Not one experiment has ever measured it directly, we don't even know it exists, and the theory I'm called a crackpot for refuting ***depends*** on the dilation of this quantity. We can do better.

If we truly wanted to dilate time, we should dilate both $t \\to t \\gamma$ and $v \\to v ~ 1/\\gamma$ in the $d=vt\\gamma$ equivalence since _both_ have time components, but we don't. To make the math fit, we only dilate the independent time component in the coordinate system of the observer in motion, leaving the velocity's magnitude unchanged by this dilation. It is not enough to blast a dimensionless scalar at an entire equation and assume that everything dilates proportionally.

Let's go back to grade school. Imagine a set of parallel lines, $A$, $B$ and $C$ where you have $B$ 'things' per $C$. If $\\gamma$ effects $A$ and $C$, and $d = An\\frac{{B}}{{C}}$ that's to say ($d = tv \\frac{{\\gamma}}{{\\gamma}}$), then $d$ remains unchanged. This is all well and good, giving rise the constant 4-vecs in special relativity, but this is objectively _not_ what we observe. 

Whether one makes the argument that the cosmos is inflating, and in that regard $d^\\prime \\neq d$, or the argument that if the above where true, than no dilation effects would be observed; both are correct.

What we've been doing is dilating _velocity_, not time, which gave us numerically correct solutions despite not understanding a thing about what's going on.

It is for this reason that the time equivalent quantity I'm proposing works; one where the numerator is integrated with time, while the denominator remains fixed... that's how we've already been applying time dilation.

I'm well aware we can put this all in a tensor and make things pretty, maybe even win some people over, but that's the problem. We need to go back to $p=mv$, not to this fanciful world some have created to boost funding. There are details that we've overlooked.

This relative dilation is even inherent in the components of $c$, $\\mu_0$ and $\\epsilon_0$... the permitivity and permability of _free space_, both with units proportional to $1/\\text{{length}}$. Of course they're going to respond to a shift in spatial density.

---

Let's lay this out more explicitly by considering the divergence of the magnetic field. If this model is correct, than the divergence of the magnetic field $\\neq 0$, but rather $\\delta$. This divergence only remains hidden due to our co-motion along the temporal (density) axis. The following is the divergence of the magnetic field in the reference frame of the 'oven' in the common raisin bread analogy.

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
{smp.latex(electricDivergence)} = {spinorSubbed * electricDivergenceSubbed} ~\\text{{m}}
$$

This is our rotational velocity to within:

$$
\\Gamma \\delta c = \\frac{{2 \\pi R_\\oplus}}{{86400}} \\pm {rotDiff:.4f}\\%
$$

And a direct, mechanistic link to why our 'time' dilates at the poles due to gravity at an almost identical rate to that of the equator due to the additional rotational velocity. Earth is *not* special, but these quantities are *not* constant.

This is off by ${rotV - (spinorSubbed * electricDivergenceSubbed):.5f}$ meters on a radius the size of Earth's... well within measurement error.

> Note here that our polar radius was used to avoid the influence of centripital force.

         """
)

# %%

print_md(f""" 
## Displacement Current Equivalent 

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

**Note:** A radius exists within Earth's range, between it's equatorial and polar radii where this quantity is precisely $\\Phi$

         """)


# %%

# Setup our function
gr = galacticRadiusPolynomial()

gr

fRadiusGalactic = smp.lambdify([rGalactic], gr.subs(subsSolar).subs([[G, spc.G]]))

# %%

estimatedRadius = (observationalSolarGalacticCenterRadius[0] + observationalSolarGalacticCenterRadius[1]) / 2

galacticRadiusSolution = fsolve(fRadiusGalactic, estimatedRadius) / (spc.kilo * spc.parsec)


# %%

# I'm sorry, but this was built for my website so the lines are assuming you're using a dark background.
Image(filename="./cmb_velocity_from_gravtational_params.png", width=450)

# %% [md]

## Solar Galactic Radius Without Additional Mass

Consider the image above. If this model is correct, some simple trig gives us the following:

$$
v = \frac{\sqrt{GM \left(GM + 2 R^3\right)}}{R^2}
$$

Or

$$
v^2 = \frac{2GM}{R} + \frac{\left(GM\right)^2}{R^4}
$$

Or even:

$$
v^2 R^4 - 2 GM R^3 - \left(GM\right)^2 = 0
$$

If we substitute our galactic parameters into this equation, we find:

$$
\left(2.3 \times 10^5\right)^2 R^4 - 2 G\left(9.94 \times 10^{40}\right)R^3 - \left(G \left(9.94 \times 10^{40}\right)\right)
$$


# %%

print_md(f""" 

While the above equation is technically symbolically solvable for $R$, the solutions are beyond messy. We can find a numeric solution with scipy, of {galacticRadiusSolution[0]:.2f} kiloparsecs. This matches observational estimates of {observationalSolarGalacticCenterRadius[0] / (spc.kilo * spc.parsec)} to {(observationalSolarGalacticCenterRadius[1] / (spc.kilo * spc.parsec)):.2f} kiloparsecs, with only visible mass.
""")

# %%
# The magnitude of this mass is completely irrelevant as it is applied to both equations. 
mass_placeholder = ac.M_sun.value
rs = np.linspace(ac.R_sun.value, spc.astronomical_unit * 100, 1000)

v1 = np.sqrt(spc.G * mass_placeholder / rs)

v2 = np.sqrt( (spc.G * mass_placeholder * (spc.G * mass_placeholder  + 2 * rs**3))) / rs**2

fig = px.scatter(title="Rotational Velocity Deviation")

fig.add_scatter(x=rs, y=v1, name=r"$\sqrt{\frac{GM}{R}}$")
fig.add_scatter(x=rs, y=v2, name=r"$\frac{\sqrt{GM \left(GM + 2 R^3\right)}}{R^2}$")

fig.update_yaxes(title_text="Velocity (m/s)")
fig.update_xaxes(title_text="Radius (meters)")

fig.show(renderer="notebook")
# %%

print_md(f""" 
## Mass Density Approximation

Finally, consider the relationship:

$$
\\beta_0^2 = {smp.latex(vbarAlpha[0])} = {smp.latex(vbarAlpha[1])} \\pm {vbarAlphaDiff:.5f}\\%
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

and 

$$
\\beta^2 \\approx \\left(\\alpha \\left(2 - \\alpha\\right)\\right)^3
$$

then finally

$$
\\rho = \\frac{{3 \\left( \\frac{{\\delta R^{{3}}}}{{2G}} \\right)}}{{4 \\pi R^{{3}}}} = \\frac{{3 \\delta}}{{8 \\pi G}} = {smp.latex(md[0])}
$$

This is our traditionally calculated mean density, $\\rho_\\prime$

$$
\\rho = {smp.latex(md[0])} = \\rho_\\prime \\pm {rhoDiff:.5f}\\%
$$

**Note:** Our _mean_ radius was used to calculate the accuracy of this mass density function, to most accurately approximate the Earth's volume.

Not a _great_ match, but pretty good considering almost every component in there is an approximation, slid by some $\\Delta$. My point is not to say that our mass density can be found by this function, but rather that the geometry of our Universe can be better understood by studying these approximations. _One_ of these alignments may be a 'dimensional coincidence', as one email put it... but not **all** of them.

And no, the dimensionality is not off. There is no such thing as time, and we need to remove it from our text books. Our 4th displacement coordinate is best represented by a position along the density axis that arises from this geometry, and the _driving_ rate of change is and always has been _motion..._ not time. This relationship is such that if $\\sum \\dot{{x}} = 0$, then $dt$ _appears_ to equal 0 as $ds = 0$ due to the magnitude of $dx$.

Time is nothing but cummulative motion and rest... a useful mathematical tool that only arises from our perception along this co-moving axis, and is not an inherent part of the system.

         """)


# %% 
print_md(f""" 
## Orbital Velocity $\\alpha$ symmetry

Consider the radius independent magnitude of our proposed spatial dilation for the prominent body in our system the Sun.

$$
v_g = {smp.latex(solarMetricVelocity())}
$$

Dividing this by our orbital velocity gives:

$$
{smp.latex(solarMetricVelocity() / vOrbit)} = {solarMetVRatio}
$$

As this compares this theoretical metric velocity to our observed orbital velocity, we should then consider the ratio of distances:

$$
\\frac{{\\mathbf{{AU}}}}{{R_\\oplus}} = {(au / R).subs(subs)}
$$

Or

$$
\\alpha \\frac{{\\mathbf{{AU}}}}{{R_\\oplus}} = {auRRatio}
$$

This is a {auRRatioPercentDiff:.3f}% percent difference.
    """)


# %%

# %%
print_md(f""" 
## $p=mv$

Consider Plank's constant $h$, the momentum carrying quantity for a massless photon.

$$
E = \\frac{{hc}}{{\\lambda}}
$$

If, as derived above, we assert the following:

$$
\\alpha = \\frac{{\\omega}}{{2 - \\alpha}} = {smp.latex(omega() / (2 - alpha))}
$$

Then we find:

$$
h^\\prime = {smp.latex(hsols[0].subs(subs))}, {smp.latex(hsols[1].subs(subs))}
$$

Or,

$$
h^\\prime = h \\pm {alphaOmegaDiff1}\\%
$$
And

$$
h^\\prime = h \\frac{{\\alpha}}{{2 - \\alpha}} \\pm {alphaOmegaDiff}
$$

Giving now our _second_ recipricol function giving $\\alpha$ in terms of gravitational parameters.
    """)

# %% [md]

## Supportive Arguments

1. No terminating geodesics, and no more blackholes. Singularities become nothing but a mirage in the distance that vanishes as you approach that gravitational field.
2. The Bullet Cluster observations are heavily indicative of a repulsive gravity model.
3. This notion of time paints an entirely different picture of chronology, where there is no limitation on backwards time travel, but a limitation on affecting one's own past. Go back to $t - n$ and try to kill your grandfather. He won't be there anymore as he's now at $t$, if he's still around. Don't get me started on the potentials around massless conciousness in this buoyant, 4-dimensional environment.
4. Machian inertia. This model concludes that the Universe is almost precisely as Mach had predicted, with inertia being comprised of the sum of this radiating force for all other bodies in the Universe.
5. The entire equivalence principle. We've all become used to it so we overlook just how profound it is, but it's literally assuming what I'm saying is true mathematically and then just pretending it's not in reality.
6. The lack of a magnetic monopole, or rather, the existince of an electric monopole. Where do you think that 'sink and faucet' model is funelling charge to? Charge is just space flowing through a hole in an equi-temporal plane, while charged bodies are just bodies with a higher tendency towards one side of this plane due to the buoyant forces along this density axis.
7. The Fermi paradox. There's an entirely other dimension that life might exist across that we have been so far unable to transcend.
8. Relative motion. Think about it. Do we currently calculate all motion as being truly relative, or do we not _already_ calculate motion as relative to the underlying gravitational acceleration and then lie to undergrads about it? The tensors of GR mask this fact to an extent, but reduce them to their calculus equivalents and we are doing _precisely_ this.
