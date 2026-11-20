from physics.symbolic.variables import c, vbar, t, gamma, tau, v, a, x, dphi, dtau, r
import sympy as smp


def variable_tau(v=v, vbar=vbar):
    return v / vbar


def dy(v=v, c=c, tau=tau):
    """Return the fifth-dimensional displacement for velocity ``v`` and ``tau``."""
    return tau * v**2 / smp.sqrt(c**2 - v**2)


def dy_dtau(v=v, c=c, tau=tau):
    """Calculate the partial derivative of ``y`` with respect to ``tau``."""
    return dy(v=v, c=c, tau=tau).diff(tau)


def dy_dx(v=v, c=c, tau=tau):
    return dy(v=v, c=c, tau=tau) / (v * tau)


def angular_momentum_L(R=r, dphi=dphi, dtau=dtau):
    return R**2 * dphi / dtau
