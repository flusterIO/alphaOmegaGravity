from importlib import import_module as _import_module

from .variables import n_r, r

from physics.symbolic.function_generators.geometry.coord_system_utils import (
    coordSysComponentFromIndex,
    radialDistanceFromOrigin,
    radialUnitVec,
    vecComponentFromIndex,
)
from physics.symbolic.function_generators.cosmological.earth import dthetaEarth
from physics.symbolic.function_generators.geometry.spheres import (
    sphereSurfaceArea,
    sphereVolume,
)
from .evaluators.evaluate_metric import check_einstein_equations
from .function_generators.cosmological.blackholes import Rs
import physics.symbolic.function_generators.electromagnetism as em
from .function_generators.electromagnetism.photon_force_ratio import (
    photon_gravity_to_electromagnetic_force_ratio,
)
import physics.symbolic.function_generators.electromagnetism.kaluza_klein as kk
import physics.symbolic.function_generators.theory.metrics as metrics
from .function_generators.theory.metrics.five_dimensional_metrics import (
    fiveDimensionalMetricTensor,
)

import physics.symbolic.function_generators.electromagnetism.four_potential as four_potential
from .function_generators.electromagnetism.four_potential import lorentz_four_force
from .function_generators.electromagnetism.four_potential.four_potential import (
    four_potential_from_electromagnetic_tensor,
)
from .function_generators.general.general_linear_algebra_quantities import vec3, vec4
from .function_generators.general.general_physical_quantities import earthMomentum
from .function_generators.geometry.dilation_matrices import (
    dilation_matrix,
    dilationMatrixRadialProjection,
    radial_dilation_metric_tensor,
)
from .function_generators.geometry.rotation_matrices import earth_rotation_matrix
from .function_generators.fluid_dynamics.fluid_equations import (
    bernoulli_pressure,
    continuity_equation,
    divergence,
    dynamic_pressure,
    gradient,
    hydrostatic_pressure,
    incompressible_continuity_equation,
    kinematic_viscosity,
    laplacian,
    mass_flow_rate,
    material_derivative,
    navier_stokes_residual,
    reynolds_number,
    vector_laplacian,
    velocity_vector,
    volumetric_flow_rate,
)
from .function_generators.geometry.vector_generatrs import vbarCosmologicalXVector
from .function_generators.newtonian.doppler import newtonianDopplerShift
from .function_generators.newtonian.gravity import (
    centripitalAcceleration,
    centripitalForce,
    dt_from_orbital_vel,
    escape_velocity,
    gravitationalPotentialEnergy,
    orbital_velocity,
    periodOfCircularOrbit,
    gravitationalAcceleration,
)
from .function_generators.quantum.alpha import alphaTraditional
from .function_generators.quantum.h import (
    plankEnergy,
    plankLength,
    plankMass,
    plankTime,
)
from .function_generators.standard_relativity.doppler import (
    gammaExpanded,
    relativisticDopplerShift,
)
from .function_generators.standard_relativity.lens import gravitationalLensDeflection
from .function_generators.standard_relativity.radial import schwartzchildRadius
import physics.symbolic.function_generators.standard_relativity.GR as GR
from .function_generators.theory.em.current_equivalent import (
    currentEquivalent,
    newtonsConstantFromAlpha,
)
from .function_generators.theory.em.electric_field_divergence import (
    electricFieldDivergence,
)
from .function_generators.theory.em.magnetic_divergence import (
    magneticDivergence,
    magneticDivergenceIntegral,
    magneticDivergenceRadialVectorMagnitude,
)
from .function_generators.theory.gravity.delta import deltaExpanded
from .function_generators.theory.gravity.gamma_equiv import (
    _Gamma,
    gammaEquivTD,
    gravitationalAccelFromElectromagnticParams,
    gammaBar,
)
from .function_generators.theory.gravity.mass_density import (
    massDensityEquivalence,
    massDensityTraditional,
    GammaFromMassDensityEquivalence,
    GammaFromRotation,
)
from physics.symbolic.function_generators.theory.gravity.bouyancy import (
    net_buoyant_force,
    object_weight,
    apparent_weight,
    buoyant_force,
    hydrostatic_pressure,
    displaced_volume_for_floating_object,
    submerged_fraction,
    floats_or_sinks_condition,
)
from .function_generators.theory.gravity.omega import omega
from .function_generators.theory.gravity.orbit import orbital_velocity_theory
from .function_generators.theory.local_parameter_expansion import (
    beta_squared_equivalent,
)
from .function_generators.theory.local_parameter_expansion.beta_squared_equivalent import (
    betaSquaredEquivalent,
)
from .function_generators.theory.quantum.alpha import alphaEquivalence, alphaSolar

from .function_generators.theory.velocity.gamma_related import gammaImaginary
from .function_generators.theory.velocity.velocity_inputs import t_r_g
from .function_generators.theory.velocity.vbar import (
    vbarAlphaEquivalence,
    vbarCosmological,
    vbarKinematic,
    vFromGamma,
)
from .operators.del_mu import del_mu
from .tensors_kinda.cristoffel import cristoffel, cristoffelTensor
from .tensors_kinda.curvature import curvatureTensor
from .tensors_kinda.F import LLF, F_sub_mu_nu, F_super_mu_nu
from .tensors_kinda.gmn import g_mn, gSuper_mn
from .variables import *
import physics.symbolic.function_generators.theory.kk.kk_theory_quantities as kk_theory
