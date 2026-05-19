import sympy as smp
import scipy.constants as spc
import astropy.units as u
import astropy.constants.iau2015 as ac
import astropy.constants.codata2010 as cd

from physics.values import unit_corrected_values


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
    au,
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
        r"\bar{h}",
        r"\l_{p}",
        r"\gamma",
        r"\mathbf{AU}",
    ]
)

_c = smp.sqrt(1 / (mo * eo))
_delta = -2 * G * M / R**3
_n = mo / (4 * smp.pi)
_omega = _delta ** smp.Rational(1, 3)
_alpha = smp.Rational(1, 2) * e**2 / h * smp.sqrt(mo / eo)
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
