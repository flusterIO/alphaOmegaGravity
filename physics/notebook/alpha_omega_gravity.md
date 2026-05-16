---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.4
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
import scipy.constants as spc
import astropy.constants.codata2010 as cd
import astropy.constants.astropyconst40 as ac
import plotly.graph_objects as go
import plotly.io as pio
import sympy as smp
import plotly.express as px
import pandas as pd
from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"
from IPython.display import display, Markdown


def print_md(string):
    display(Markdown(string))
```

## Parameters

As I'm offline while writing this I'm unable to look up the papers these parameters were taken from, but all values were taken from semi-recent, well cited sources.

```python
cmbRelativeVelocities = [
    1.228 * 10**-3 * cd.c,
    1.231 * 10**-3 * cd.c,
    1.234 * 10**-3 * cd.c,
]


earthRadiusPolar: float = 6356752

earthRadiusMean: float = 6371000
```

```python
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
]

subsSolar = [[R, ac.R_sun.base.tolist()], [M, ac.M_sun.base.tolist()]]

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
]
```

```python
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
```

```python
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
vbarKinematicDiff = percentDiff(cmbRelativeVelocities[2], vbarKinematicSubbed * u.m / u.s)
vbarKinematicDiff1 = percentDiff(cmbRelativeVelocities[0], vbarKinematicSubbed * u.m / u.s)
rhoDiff = percentDiff(ac.M_earth.value.tolist() / (4/3 * spc.pi * earthRadiusMean**3), rhoSubbed)
jDiff = abs(percentDiff(_delta.subs(subs) * _vbarCosmological.subs(subs), spc.golden_ratio))
currentEquivDiff = abs(percentDiff(currentEquivalent().subs(subs), spc.golden_ratio))
alphaSolarDiff = percentDiff(spc.alpha + 1, _alphaSolar[0].subs(subsSolar).subs(subs))
rotDiff = percentDiff(rotV, spinorSubbed * electricDivergenceSubbed)
omegaAlphaDiff = percentDiff(alphaFactor.subs([[alpha, spc.alpha]]), omegaSubbed.evalf())
```

```python
df = pd.DataFrame([ 
    {
        "derivation": "cosmological velocity from local gravitational parameters",
        "function": "vbarCosmological",
        "% error": "Difficult to assess"
    },
    {
        "derivation": "kinematic velocity from local gravitational parameters",
        "function": "vbarKinematic",
        "% error": vbarKinematicDiff
    },
    {
        "derivation": "Earth's roational velocity via $\\alpha$",
        "function": "rotationVelocityFromDivergence",
        "% error": rotDiff
    },
    {
        "derivation": "α (solar)",
        "function": "alphaSolar",
        "% error": alphaSolarDiff
    },
    {
        "derivation": "α (Earth)",
        "function": "omega, alphaFactor",
        "% error": omegaAlphaDiff
    },
    {
        "derivation": "α (from derived velocity)",
        "function": "vbarAlphaEquivalence",
        "% error": vbarAlphaDiff
    },
    {
        "derivation": "J ≈ Φ",
        "function": "currentEquivalent",
        "% error": currentEquivDiff
    },
    {
        "derivation": "Φ mass density approximation",
        "function": "massDensityEquivalence",
        "% error": rhoDiff
    },
])

df
```

<!-- #region -->
## Divergent Gravity

Here is a short derivation, but for a longer derivation, please see the documentation at [flusterapp.com](https://flusterapp.com/docs). There is a more complete derivation included as an example for the application.


### The Lighthouse & The Clocktower

> A quick thought experiment to justify a case for classical synchronicity

Consider two observers, $A$ and $B$. Let $A$ and $B$ agree that $B$ should travel between two arbitrary points in space (assuming a uniform gravitational field), $p_1$ and $p_2$ at an _agreed_ upon velocity $v$. 

At $p_2$ place a time keeping device that remains visible from $p_1$. Allow that $A$ remains at rest relative to $p_1$, and very near $p_1$ while $B$ goes in motion (at an inertial pace) between $p_1$ and $p_2$. 

At what time does $A$ observe $B$ reach $p_2$? If the value is not $d/v$, but $\gamma d/v$, then $v \neq v_\prime$. 

The only possible solution that allows for Einstein's additional factor of $\gamma$ while allowing $A$ to remain in an unchanged state is that $d \to d\gamma$ in the coordinate system of the observer in relative motion, while $v \to v/\gamma$ in the coordinate system of the observer at relative rest to accomodate the change in spatial density. Any other solution would modify the agreed upon belocity, or result in a disagreement on the time that $B$ would reach $p_2$.

If the two observers are to agree upon the time at which $B$ reaches $p_2$ while $v = v_\prime$ in the coordinate system of $A$ (the original coordinate system) then distance _must_ elongate with time.

Think about this even further. Do we not already believe this? We believe cosmic inflation is a time dependent property, what is the big leap to considering that maybe this property is _velocity_ dependent? Just look at the way the physics aligns when we make this single modification, while every single experimental validation of either SR or GR is maintained due to the $d=vt\gamma$ equivalence that both models seek to maintain (and their dependence on the equivalence principle).

It's time to think for yourself people... Einstein's not here to hold our hand anymore.
<!-- #endregion -->

```python
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
```

```python
print_md(f""" 
## Electromagnetism

If this model is correct, and Earth dilates _with_ the rest of the Universe, along with all other bodies, than we should consider the divergence of $\\delta$... or rather, $\\delta$ **as** divergence.

Taking the inverse of $\\nabla \\cdot \\omega = \\delta$, assuming a diagonal 3-dimensional matrix gives:

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

$\\alpha$ is _not_ a constant in itself, and it needs to stop being viewed as such. It's comprised of other constants, each with it's own properties and antagonists, giving $\\alpha$ different units in different contexts.

$$
\\overbrace{{\\alpha}}^{{dt}} = \\underbrace{{\\frac{{1}}{{2 - \\alpha}}}}_{{1/dx}} \\overbrace{{\\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}}}}^{{ds}}
$$

Consider the function for our Sun, where our 'Spinor' approximation limits at:

$$
\\lim_{{\\Gamma \\to (1/2)^{{1/3}}}} \\frac{{1}}{{2 - \\alpha}} \\sqrt[3]{{2\\frac{{GM_\\oplus}}{{R_\\oplus^3}}}} \\approx 1 \\pm \\alpha
$$

Or, from Python and Astropy:

$$
\\alpha_{{\\text{{Solar}}}} = {smp.latex(_alphaSolar)} = {_alphaSolar[0].subs(subsSolar).subs(subs):.8f}, {_alphaSolar[1].subs(subsSolar).subs(subs):.8f} = 1 \\pm \\alpha \\pm {alphaSolarDiff:.4f}\\%
$$

> Ignore the goofy exponents if they appear astronomical. That's coming directly from sympy, but they evaluate to 1/3.

Above, 'dt' is used simply to explain the phenomena in terms we are familiar with. There is absolutely no need for a temporal axis, and it plays absolutely no role in the underlying physics. It is simply a dimensionless ratio of change in the form of $ds/dx$, giving rise to the **sensation** of time due to our comotion and the dilation of $c$ along this axis. There is no drag, no resistance or boost... no physical interaction with the 'time medium' of any kind. There is no such thing.
""")
```

```python
print_md(
    f"""
Riddle me this... how can $c$ be constant in-between reference frames, when it's a time dependent quantity and time dilates? This logic is self defeating. _If_ $c$ is of an equivalent magnitude ***in*** each reference frame, than it ***must*** dilate between reference frames to accomodate the shift in unit magnitude. It cannot possibly be both, and this is why $c$ _must_ dilate between reference frames. It's inherent in the components of $c$, $\\mu_0$ and $\\epsilon_0$... the permitivity and permability of _free space_, both with units proportional to $1/\\text{{length}}$. Of course they're going to respond to a shift in spatial density.


Let's lay this out more explicitly, since many academics remain determined to pray to their volcano.

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
{smp.latex(electricDivergence)} = {spinorSubbed * electricDivergenceSubbed} ~\\text{{m}}
$$

This is our rotational velocity to within:

$$
\\frac{{2 \\pi R_\\oplus}}{{86400}} \\pm {rotDiff:.4f}\\%
$$

And a direct, mechanistic link to why our 'time' dilates at the poles due to gravity at an almost identical rate to that of the equator due to the additional rotational velocity. Earth is *not* special, but these quantities are *not* constant.

This is off by ${rotV - (spinorSubbed * electricDivergenceSubbed):.5f}$ meters... well within measurement error.

> Note here that our polar radius was used to avoid the influence of centripital force.

         """
)
```

```python
print_md(f""" 
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

         """)
```

```python
print_md(f""" 
Finally, consider the relationship:

$$
v_0^2 = {smp.latex(vbarAlpha[0])} = {smp.latex(vbarAlpha[1])} \\pm {vbarAlphaDiff:.5f}\\%
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

This is our traditionally calculated mean density, $\\rho_\\prime$

$$
\\rho_\\prime \\pm {rhoDiff:.5f}\\%
$$

Not a _great_ match, but pretty good considering almost every component in there is an approximation, slid by some $\\Delta$. My point is not to say that our mass density can be found by this function, but rather that the geometry of our Universe can be better understood by studying these approximations. _One_ of these alignments may be a 'dimensional coincidence', as one email put it... but not **all** of them.

And no, the dimensionality is not off. There is no such thing as time, and we need to remove it from our text books. Our 4th displacement coordinate is best represented by a position along the density axis that arises from this geometry, and the _driving_ rate of change is and always has been _motion..._ not time.

Time is nothing but cummulative motion and rest... a useful mathematical tool that only arises from our perception along this co-moving axis, and is not an inherent part of the system.

         """)
```
