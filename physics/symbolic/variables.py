from dataclasses import dataclass

import astropy.constants.codata2010 as cd
import astropy.constants.iau2015 as ac
import astropy.units as u
import numpy as np
import scipy.constants as spc
import sympy as smp
from sympy.physics.vector import ReferenceFrame
from sympy.vector import CoordSys3D, matrix_to_vector
import physics.data.li_observation as li

from physics.values import unit_corrected_values

A = CoordSys3D("A")
S = ReferenceFrame("S")
B = ReferenceFrame("B")
E = ReferenceFrame("E")


Ax, Ay, Az = A.base_scalars()
Ex, Ey, Ez = E[0], E[1], E[2]
Bx, By, Bz = B[0], B[1], B[2]
Sx, Sy, Sz = S[0], S[1], S[2]

phi = smp.Function(
    r"\phi",
)(Sx, Sy, Sz)

# Radial coordinate and refractive-index profile used by geometric-optics models.
r = smp.symbols("r", positive=True)
R = smp.symbols("R", positive=True)
n_r = smp.Function("n")(r)
_path_parameter = smp.symbols("_path_parameter", real=True)


t_r = smp.Function("t")(r)


(
    x,
    y,
    z,
    T0,
    U,
    G,
    q,
    qtest,
    g,
    kappa,
    depth,
    rho_object,
    theta,
    thetaInclination,
    M,
    Mearth,
    m,
    Rcal,
    Rearth,
    Rhat,
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
    dphi,
    dtau,
    dt,
    alpha,
    h,
    height,
    hbar,
    lp,
    gamma,
    gammaBar,
    Gamma,
    v,
    vx,
    vy,
    vz,
    V,
    vbar,
    beta,
    b,
    TDay,
    t,
    AU,
    rho,
    varrho,
    ux,
    uy,
    uz,
    ut,
    Dcal,
    Lcal,
    L,
    vgroup,
    mu,
    nu,
    Phi,
    varphi,
    S4,
    S0,
    tau,
    dy,
    dx,
    a,
    xdot,
    xhat,
) = smp.symbols(
    [
        "x",
        "y",
        "z",
        r"\mathcal{T}",
        "U",
        "G",
        "q",
        r"q_{\text{test}}",
        "g",
        r"\kappa",
        r"d_{\text{depth}}",
        r"\rho_{\text{object}}",
        r"\theta",
        r"\theta_{\scriptsize\text{inclination}}",
        "M",
        r"M_{\oplus}",
        "m",
        r"\mathcal{R}",
        r"R_{\oplus}",
        r"\hat{R}",
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
        r"\mathrm{d}\theta",
        r"\mathrm{d}\phi",
        r"\mathrm{d}\tau",
        r"\mathrm{d}t",
        r"\alpha",
        "h",
        r"h_{\text{height}}",
        r"\bar{h}",
        r"\l_{p}",
        r"\gamma",
        r"\bar{\gamma}",
        r"\Gamma",
        "v",
        r"v_x",
        r"v_y",
        r"v_z",
        "V",
        r"\bar{v}",
        r"\beta",
        "b",
        r"T_{\text{day}}",
        "t",
        r"\mathbf{AU}",
        r"\rho",
        r"\varrho",
        "u_{x}",
        "u_{y}",
        "u_{z}",
        "u_{t}",
        r"\mathcal{D}",
        r"\mathcal{L}",
        "L",
        r"v_{\text{group}}",
        r"\mu",
        r"\nu",
        r"\Phi",
        r"\varphi",
        r"S_{4}",
        r"S_{0}",
        r"\tau",
        r"\mathrm{d}y",
        r"\mathrm{d}x",
        "a",
        r"\dot{x}",
        r"\widehat{x}",
    ]
)


vOfGamma = smp.Function(
    "v",
)(gamma)

p = smp.Function("p")(x, y, z, t)
rho_field = smp.Function(r"\rho")(x, y, z, t)
u_x = smp.Function(r"u_x")(x, y, z, t)
u_y = smp.Function(r"u_y")(x, y, z, t)
u_z = smp.Function(r"u_z")(x, y, z, t)
f_x = smp.Function(r"f_x")(x, y, z, t)
f_y = smp.Function(r"f_y")(x, y, z, t)
f_z = smp.Function(r"f_z")(x, y, z, t)

_c = smp.sqrt(1 / (mo * eo))
_delta = -2 * G * M / R**3
_n = mo / (4 * smp.pi)
_omega = (-_delta) ** smp.Rational(1, 3)
_alpha = smp.Rational(1, 2) * e**2 / h * smp.sqrt(mo / eo)
_mass_plank = smp.sqrt(h * c / (2 * smp.pi * G))
_lp = smp.sqrt(hbar * G / c**3)

subs = [
    [c, spc.c],
    [G, spc.G],
    [Rearth, ac.R_earth.base.tolist()],
    [Mearth, ac.M_earth.base.tolist()],
    [me, spc.m_e],
    [mo, cd.mu0.base.tolist()],
    [eo, cd.eps0.base.tolist()],
    [e, cd.e.base.tolist()],
    [h, cd.h.value.tolist()],
    [AU, spc.au],
]


subsEarth = [
    [thetaInclination, smp.rad(23.439281)],
    [M, float(ac.M_earth.value)],
    [R, float(ac.R_earth.value)],
    [Gamma, 1 / (2 - alpha)],
    [dtheta / dt, 2 * smp.pi / TDay],
    [TDay, spc.day],
]

subsSolar = [[R, ac.R_sun.base.tolist()], [M, ac.M_sun.base.tolist()]]

subsLi = [[L, 20_000], [vgroup, li.GROUP_VELOCITY]]

expand_g = [[g, G * M / R**2]]
expand_gammaBar = [[gamma, 1 + G * M / R**3]]
expand_vbar = [[vbar, c * smp.sqrt(1 - 1 / gamma**2)]]

expand = [
    [gammaBar, 1 + g / R],
    [vOfGamma, c * smp.sqrt(1 - 1 / gamma**2)],
    *expand_g,
    *expand_vbar,
    [n, _n],
    [c, _c],
    [omega, _omega],
    [delta, _delta],
    [alpha, _alpha],
    [m_plank, _mass_plank],
    [lp, _lp],
    [hbar, h / (smp.pi * 2)],
]

expandRho = [[rho, M / (smp.Rational(4, 3) * smp.pi * R**3)]]


subsAll = [*subs, *subsEarth, *expand, *subs]


@dataclass
class Subs:
    subs = subsAll
    subsSolar = subsSolar
    expand = expand

    def subsWithoutVariable(self, var):
        """
        Uses subsAll subs.
        """
        items = []
        for v in self.subs:
            if v[0] != var:
                items.append(v)
        return items

    def expandWithoutVariable(self, var):
        items = []
        for v in self.expand:
            if v[0] != var:
                items.append(v)
        return items

    def subsSolarWithoutVariable(self, var):
        items = []
        for v in self.subsSolar:
            if v[0] != var:
                items.append(v)
        return items

    def subRHat(self, coords: ReferenceFrame, radial_unit_vec=None):
        x = coords[0]
        y = coords[1]
        z = coords[2]
        r = smp.sqrt(x**2 + y**2 + z**2)
        return [[Rhat, (x * coords.x + y * coords.y + z * coords.z) / r]]

    def subRadiusPoints(
        self,
        coords: (
            ReferenceFrame | tuple[smp.Symbol, smp.Symbol, smp.Symbol] | None
        ) = None,
        r=1,
        radial_unit_vec=None,
        radius_variable=R,
    ):
        if radial_unit_vec is None:
            radial_unit_vec = np.array([1, 0, 0])
        if coords is None:
            coords = (x, y, z)
        return [
            [coords[0], radial_unit_vec[0] * r],
            [coords[1], radial_unit_vec[1] * r],
            [coords[2], radial_unit_vec[2] * r],
        ]

    @staticmethod
    def radialVecFromPhiTheta(theta: float, phi: float):
        return np.array(
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
        )

    def subRadius(self, coords: ReferenceFrame, radial_unit_vec=None):
        x = coords[0]
        y = coords[1]
        z = coords[2]
        if radial_unit_vec is None:
            radial_unit_vec = np.array([1, 0, 0])
        return [
            [R, smp.sqrt(x**2 + y**2 + z**2)],
            [Rhat, smp.Matrix(Rearth * radial_unit_vec)],
        ]
