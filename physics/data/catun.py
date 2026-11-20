"""
Baryonic mass model of the Milky Way from Cautun et al. (2020), MNRAS 494, 4291 (arXiv:1911.04557).

Every number below is transcribed from the paper:
  Table 1  (fixed parameters)          -> BULGE_FIXED, DISC_FIXED, GAS_FIXED, CGM
  Table 2  (best-fit parameters)       -> SETS['contracted'] (their preferred model) and SETS['nfw']
  Eq. (1)  bulge, Eq. (3) stellar discs, Eq. (4) HI/H2 discs, Eq. (5)-(6) CGM

Units: kpc, Msun, km/s.  G = 4.30091e-6 kpc (km/s)^2 / Msun.
Two kinds of output:
  * enclosed_mass(r)   : spherical enclosed mass M_b(<r) of every component (what the paper's
                         v_orbit = sqrt(2 G M/R) uses)
  * v2_newton(R)       : in-plane Newtonian circular speed^2 of the baryons, with the flattened discs
                         treated properly (layered-ring potential). Used as g_bar for the RAR-type test.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import integrate, special

G = 4.30091e-6  # kpc (km/s)^2 / Msun
FBAR = 0.157  # cosmic baryon fraction used by Cautun et al. (their Sec. 2.4)

# ---------------- Table 1: fixed parameters ----------------
BULGE_FIXED = dict(r0=0.075, rcut=2.1, alpha=1.8, q=0.5)  # kpc
DISC_FIXED = dict(thin=dict(zd=0.300), thick=dict(zd=0.900))  # kpc
GAS_FIXED = dict(
    HI=dict(zd=0.085, Rm=4.0, Rd=7.0, S0=53.0),  # S0 in Msun/pc^2
    H2=dict(zd=0.045, Rm=12.0, Rd=1.5, S0=2200.0),
)
CGM = dict(A=0.190, beta=-1.46)

# ---------------- Table 2: best-fit parameters ----------------
SETS = {
    "contracted": dict(
        rho0_bulge=103.0,  # Msun/pc^3
        S0_thin=731.0,
        Rd_thin=2.63,  # Msun/pc^2, kpc
        S0_thick=101.0,
        Rd_thick=3.80,
        M200_tot=1.08e12,
        R200=218.0,
    ),  # Msun, kpc   (for the CGM)
    "nfw": dict(
        rho0_bulge=101.0,
        S0_thin=1070.0,
        Rd_thin=2.43,
        S0_thick=113.0,
        Rd_thick=3.88,
        M200_tot=0.95e12,
        R200=207.0,
    ),
}


# ---------------- surface-density helpers ----------------
def sigma_exp(a, S0, Rd):  # stellar discs, Msun/kpc^2
    return S0 * 1e6 * np.exp(-a / Rd)


def sigma_hole(a, S0, Rm, Rd):  # gas discs (with inner hole), Msun/kpc^2
    return S0 * 1e6 * np.exp(-Rm / a - a / Rd)


def vertical_profile(kind, zd):
    """h(z) normalised so that int h dz = 1 over (-inf,inf)."""
    if kind == "exp":
        return lambda z: np.exp(-np.abs(z) / zd) / (2 * zd)  # Eq. (3)
    if kind == "sech2":
        return lambda z: 1.0 / np.cosh(z / (2 * zd)) ** 2 / (4 * zd)  # Eq. (4)


def components(p, stellar_scale=1.0, cgm_scale=1.0):
    """Return list of dicts describing every baryonic component of parameter set p."""
    c = []
    c.append(
        dict(
            name="bulge",
            type="bulge",
            rho0=p["rho0_bulge"] * 1e9 * stellar_scale,
            **BULGE_FIXED,
        )
    )
    c.append(
        dict(
            name="thin",
            type="disc",
            kind="exp",
            zd=DISC_FIXED["thin"]["zd"],
            sigma=lambda a, p=p, s=stellar_scale: s
            * sigma_exp(a, p["S0_thin"], p["Rd_thin"]),
        )
    )
    c.append(
        dict(
            name="thick",
            type="disc",
            kind="exp",
            zd=DISC_FIXED["thick"]["zd"],
            sigma=lambda a, p=p, s=stellar_scale: s
            * sigma_exp(a, p["S0_thick"], p["Rd_thick"]),
        )
    )
    for nm in ("HI", "H2"):
        g = GAS_FIXED[nm]
        c.append(
            dict(
                name=nm,
                type="disc",
                kind="sech2",
                zd=g["zd"],
                sigma=lambda a, g=g: sigma_hole(a, g["S0"], g["Rm"], g["Rd"]),
            )
        )
    c.append(
        dict(name="CGM", type="cgm", M200=p["M200_tot"] * cgm_scale, R200=p["R200"])
    )
    return c


# ---------------- spherical enclosed masses ----------------
def bulge_enclosed(r, comp):
    rho = (
        lambda x: comp["rho0"]
        / (1 + x / comp["r0"]) ** comp["alpha"]
        * np.exp(-((x / comp["rcut"]) ** 2))
    )
    # ellipsoid r'<r (contained in the sphere of radius r because q<1): dV = 4 pi q r'^2 dr'
    return (
        4
        * np.pi
        * comp["q"]
        * integrate.quad(
            lambda x: x**2 * rho(x), 0, r, limit=200, points=[comp["r0"], comp["rcut"]]
        )[0]
    )


def disc_enclosed(r, comp):
    h = vertical_profile(comp["kind"], comp["zd"])
    sig = comp["sigma"]
    f = lambda R, z: 2 * np.pi * R * sig(R) * h(z)
    val = integrate.dblquad(
        f,
        0,
        r,
        lambda z: 0.0,
        lambda z: np.sqrt(max(r * r - z * z, 0.0)),
        epsabs=1e3,
        epsrel=1e-6,
    )[0]
    return 2 * val  # z>0 and z<0


def cgm_enclosed(r, comp):  # Eq. (6)
    A, b = CGM["A"], CGM["beta"]
    return 3 * A / (b + 3) * FBAR * comp["M200"] * (r / comp["R200"]) ** (b + 3)


def enclosed_mass(r, p, **kw):
    out = {}
    for c in components(p, **kw):
        out[c["name"]] = {
            "bulge": bulge_enclosed,
            "disc": disc_enclosed,
            "cgm": cgm_enclosed,
        }[c["type"]](r, c)
    out["total"] = sum(out.values())
    return out


# ---------------- in-plane Newtonian circular speed (flattened discs) ----------------
_zn, _zw = np.polynomial.legendre.leggauss(24)


def _layer_potential(R, comp):
    """Potential at (R, z=0) of a disc with surface density sigma(a) and vertical profile h(z)."""
    sig = comp["sigma"]
    h = vertical_profile(comp["kind"], comp["zd"])
    zd = comp["zd"]
    zmax = 12 * zd
    zs = 0.5 * zmax * (_zn + 1)
    ws = 0.5 * zmax * _zw  # z in (0,zmax)
    tot = 0.0
    for z, w in zip(zs, ws):

        def integrand(a):
            d2 = (a + R) ** 2 + z * z
            m = 4 * a * R / d2
            return (
                2
                * np.pi
                * a
                * sig(a)
                * special.ellipk(min(m, 1 - 1e-14))
                * 2
                / (np.pi * np.sqrt(d2))
            )

        val = integrate.quad(
            integrand, 1e-4, 80.0, points=[R], limit=300, epsabs=1e-3, epsrel=1e-8
        )[0]
        tot += 2 * w * h(z) * val  # 2: both z signs
    return -G * tot


def v2_newton(R, p, **kw):
    """Baryonic v_c^2 (km/s)^2 at radius R [kpc], bulge+CGM treated as spherical, discs as layered rings."""
    comps = components(p, **kw)
    v2 = 0.0
    parts = {}
    for c in comps:
        if c["type"] == "disc":
            d = 1e-3 * R
            gR = (
                -(_layer_potential(R + d, c) - _layer_potential(R - d, c))
                / (2 * d)
                * (-1)
            )  # g = -dPhi/dR
            gR = (_layer_potential(R + d, c) - _layer_potential(R - d, c)) / (
                2 * d
            )  # dPhi/dR >0 -> inward pull
            parts[c["name"]] = R * gR
        else:
            m = {"bulge": bulge_enclosed, "cgm": cgm_enclosed}[c["type"]](R, c)
            parts[c["name"]] = G * m / R
    parts["total"] = sum(parts.values())
    return parts


def getDataFrame():
    return pd.read_csv(Path(__file__).parent / "tabular" / "catun_observation.csv")


if __name__ == "__main__":
    p = SETS["contracted"]
    for r in (8.12,):
        m = enclosed_mass(r, p)
        print(
            f"M_b(<{r} kpc) by component [1e10 Msun]:",
            {k: round(v / 1e10, 3) for k, v in m.items()},
        )
    tot = {c["name"]: None for c in components(p)}
    print(
        "component total masses [1e10 Msun] (thin disc analytic 2 pi S0 Rd^2 = %.2f)"
        % (2 * np.pi * p["S0_thin"] * 1e6 * p["Rd_thin"] ** 2 / 1e10)
    )
    print(
        "bulge total:",
        round(bulge_enclosed(40.0, components(p)[0]) / 1e10, 3),
        " (paper: 0.94)",
    )
