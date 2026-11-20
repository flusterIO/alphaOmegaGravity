from physics.symbolic.variables import ux, uy, uz, ut
import sympy as smp


def vec3(_ux=ux, _uy=uy, _uz=uz):
    smp.Matrix([_ux, _uy, _uz])


def vec4(_ut=ut, _ux=ux, _uy=uy, _uz=uz):
    smp.Matrix([_ut, _ux, _uy, _uz])
