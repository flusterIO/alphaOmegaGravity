"""Published bulk-flow and observer-velocity measurements.

Quality scores are editorial judgments (1 = weak or highly disputed evidence;
10 = exceptionally precise, well-characterized measurement). They are not
scores assigned by the cited authors. Scores consider sample size, statistical
and systematic uncertainty, analysis transparency, and independent validation.
Measurements sharing source data are identified in ``notes`` and should not be
treated as independent.
"""

from enum import Enum

import pandas as pd


class VelocityType(Enum):
    """Physical quantity represented by a reported speed."""

    COSMOLOGICAL = "cosmological bulk or peculiar flow"
    KINEMATIC = "observer motion relative to the CMB"


class VelocityMethod(Enum):
    """Primary data or comparison used to infer the velocity."""

    SNE = "Type Ia/II supernova distances and redshifts"
    SNE_CORRELATED_PECULIAR_VELOCITIES = "Correlated Type Ia supernova peculiar-velocity covariance"
    CMB_DIPOLE = "CMB temperature dipole"
    CMB_DOPPLER_ABERRATION = "CMB Doppler boosting and aberration"
    CMB_INTRINSIC_DIPOLE_TEST = "CMB Doppler/aberration, allowing intrinsic dipole"
    CMB_KSZ = "CMB kinetic Sunyaev-Zeldovich signal from clusters/galaxies"
    SNE_CMB_COMPARISON = "Supernova flow constrained or compared with CMB direction/frame"
    SNE_PREDICTED_FLOW_COMPARISON = "Supernova velocities compared with galaxy-density flow model"
    MULTI_TRACER = "Combined peculiar-velocity distance-indicator catalogs"
    CLUSTER_DISTANCES = "Galaxy-cluster distance indicators"


class ResultKind(Enum):
    """How to interpret the reported value and uncertainty fields."""

    ESTIMATE = "estimate"
    UPPER_LIMIT = "upper limit"
    REPORTED_RANGE = "reported range"
    COSMOLOGICAL_PARAMETER = "cosmological parameter constraint, not a speed"


_MEASUREMENTS = (
    {
        "paper": "Cosmological constraints from Type Ia supernovae peculiar velocity measurements",
        "source_url": "https://doi.org/10.1103/PhysRevLett.99.081301",
        "velocity_km_s": None,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE_CORRELATED_PECULIAR_VELOCITIES,
        "result_kind": ResultKind.COSMOLOGICAL_PARAMETER,
        "margin_of_error_km_s": None,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby Type Ia SNe",
        "effective_scale_mpc_h": None,
        "year": 2007,
        "quality_score": 8,
        "notes": "Gordon, Land & Slosar detect correlated SNe peculiar velocities and infer sigma_8 = 0.79 ± 0.22; the paper does not report a single bulk-flow speed.",
    },
    {
        "paper": "Determining the motion of the Solar system relative to the CMB using Type Ia supernovae (uncorrelated-PV fit)",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2008.13239.x",
        "velocity_km_s": 475,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.SNE_CMB_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 134,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby Type Ia SNe",
        "effective_scale_mpc_h": None,
        "year": 2008,
        "quality_score": 5,
        "notes": "Gordon, Land & Slosar estimate Solar System motion relative to the CMB while omitting correlated peculiar-velocity covariance. The 525.49 km/s comparison value lies within the reported 1σ interval; authors find this fit overstates detection significance.",
    },
    {
        "paper": "Determining the motion of the Solar system relative to the CMB using Type Ia supernovae (correlated-PV fit)",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2008.13239.x",
        "velocity_km_s": 468,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.SNE_CMB_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 186,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby Type Ia SNe",
        "effective_scale_mpc_h": None,
        "year": 2008,
        "quality_score": 7,
        "notes": "Preferred covariance-aware estimate of Solar System motion relative to the CMB. The 525.49 km/s comparison value is inside the reported 1σ interval; including correlations reduces the detection from about 3.5σ to 2.5σ.",
    },
    {
        "paper": "Determining the motion of the Solar system relative to the CMB using Type Ia supernovae (uncorrelated-PV Local Group fit)",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2008.13239.x",
        "velocity_km_s": 697,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.SNE_CMB_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 137,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby Type Ia SNe",
        "effective_scale_mpc_h": None,
        "year": 2008,
        "quality_score": 5,
        "notes": "Local Group motion inferred by Gordon, Land & Slosar after subtracting the Solar System motion relative to the Local Group; uncorrelated-PV analysis.",
    },
    {
        "paper": "Determining the motion of the Solar system relative to the CMB using Type Ia supernovae (correlated-PV Local Group fit)",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2008.13239.x",
        "velocity_km_s": 690,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.SNE_CMB_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 201,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby Type Ia SNe",
        "effective_scale_mpc_h": None,
        "year": 2008,
        "quality_score": 7,
        "notes": "Local Group motion inferred from the covariance-aware Solar System estimate by subtracting Solar System–Local Group motion; not the same quantity as the Solar System speed alone.",
    },
    {
        "paper": "Measuring the cosmological bulk flow using the peculiar velocities of supernovae",
        "source_url": "https://arxiv.org/abs/1102.0800",
        "velocity_km_s": 188,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": (-103, 119),
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "z < 0.05",
        "effective_scale_mpc_h": None,
        "year": 2011,
        "quality_score": 7,
        "notes": "68% interval; SNe Ia bulk-flow fit. The authors reject the then-claimed >600 km/s dark flow at >99.7% confidence.",
    },
    {
        "paper": "Cosmic flows in the nearby universe from Type Ia supernovae",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2011.20050.x",
        "velocity_km_s": 249,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE_PREDICTED_FLOW_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 76,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": "50 (minimum-variance Gaussian window)",
        "year": 2012,
        "quality_score": 7,
        "notes": "245 SNe Ia; also compares the measured flow with PSCz galaxy-density predictions. The SNe and model-residual results are listed separately.",
    },
    {
        "paper": "Cosmic flows in the nearby universe from Type Ia supernovae (PSCz residual)",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2011.20050.x",
        "velocity_km_s": 150,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE_PREDICTED_FLOW_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 43,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": None,
        "year": 2012,
        "quality_score": 7,
        "notes": "Residual SNe bulk motion not accounted for by the PSCz density-field prediction; not a direct CMB-dipole measurement.",
    },
    {
        "paper": "Measuring cosmic bulk flows with Type Ia supernovae from the Nearby Supernova Factory (0.045 < z < 0.06)",
        "source_url": "https://doi.org/10.1051/0004-6361/201321880",
        "velocity_km_s": 650,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 398,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "0.045 < z < 0.06",
        "effective_scale_mpc_h": None,
        "year": 2014,
        "quality_score": 5,
        "notes": "Low-significance shell estimate with a large uncertainty; the paper reports no significant bulk flow in its higher-redshift shell.",
    },
    {
        "paper": "Measuring cosmic bulk flows with Type Ia supernovae from the Nearby Supernova Factory (CMB-direction, high-z limit)",
        "source_url": "https://doi.org/10.1051/0004-6361/201321880",
        "velocity_km_s": None,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE_CMB_COMPARISON,
        "result_kind": ResultKind.UPPER_LIMIT,
        "margin_of_error_km_s": None,
        "upper_limit_km_s": 240,
        "reported_range_km_s": None,
        "redshift_range": "0.06 < z < 0.1",
        "effective_scale_mpc_h": None,
        "year": 2014,
        "quality_score": 7,
        "notes": "68% upper limit with the flow direction fixed to the CMB dipole; this is not a point measurement of 240 km/s.",
    },
    {
        "paper": "Planck 2018 results. I. Overview and the cosmological legacy of Planck",
        "source_url": "https://doi.org/10.1051/0004-6361/201833910",
        "velocity_km_s": 369.82,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.CMB_DIPOLE,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 0.11,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": None,
        "year": 2020,
        "quality_score": 10,
        "notes": "Solar System barycentre relative to the CMB; the precision assumes the observed dipole is kinematic.",
    },
    {
        "paper": "Planck 2013 results. XXVII. Doppler boosting of the CMB: Eppur si muove",
        "source_url": "https://doi.org/10.1051/0004-6361/201321556",
        "velocity_km_s": 384,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.CMB_DOPPLER_ABERRATION,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": {"statistical": 78, "systematic": 115},
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": None,
        "year": 2014,
        "quality_score": 7,
        "notes": "Independent Doppler/aberration signature estimate; statistical and systematic errors are kept separate.",
    },
    {
        "paper": "First constraints on the intrinsic CMB dipole and our velocity with Doppler and aberration",
        "source_url": "https://doi.org/10.1103/PhysRevLett.127.101301",
        "velocity_km_s": 300,
        "velocity_type": VelocityType.KINEMATIC,
        "method": VelocityMethod.CMB_INTRINSIC_DIPOLE_TEST,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": (-93, 111),
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": None,
        "year": 2021,
        "quality_score": 8,
        "notes": "SMICA result; allows a possible intrinsic dipole rather than assigning the full CMB temperature dipole to motion.",
    },
    {
        "paper": "First measurement of the bulk flow of nearby galaxies using the cosmic microwave background",
        "source_url": "https://doi.org/10.1093/mnras/stt525",
        "velocity_km_s": 533,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.CMB_KSZ,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 263,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": "50",
        "year": 2013,
        "quality_score": 6,
        "notes": "Kinetic Sunyaev-Zeldovich estimate for nearby galaxies; large uncertainty and an independent method from distance-indicator flows.",
    },
    {
        "paper": "Consistently large cosmic flows on scales of 100 h^-1 Mpc",
        "source_url": "https://doi.org/10.1111/j.1365-2966.2008.14089.x",
        "velocity_km_s": 407,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.MULTI_TRACER,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 81,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": "50 (Gaussian window)",
        "year": 2009,
        "quality_score": 7,
        "notes": "Minimum-variance combination of peculiar-velocity surveys. The authors discuss sparse sampling and window-function effects.",
    },
    {
        "paper": "Cosmic flows on 100 Mpc/h scales: standardized minimum variance bulk flow, shear and octupole moments",
        "source_url": "https://arxiv.org/abs/0911.5516",
        "velocity_km_s": 416,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.MULTI_TRACER,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 78,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": None,
        "effective_scale_mpc_h": "100",
        "year": 2010,
        "quality_score": 7,
        "notes": "Minimum-variance reanalysis of largely overlapping peculiar-velocity catalogs with shear/octupole control; not independent of Watkins et al. (2009).",
    },
    {
        "paper": "Peculiar-velocity cosmology with Types Ia and II supernovae",
        "source_url": "https://doi.org/10.1093/mnras/stab1446",
        "velocity_km_s": 245,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": (-31, 32),
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "low-z supernova sample",
        "effective_scale_mpc_h": "~30",
        "year": 2021,
        "quality_score": 8,
        "notes": "Bulk flow from a 775-object compilation of low-redshift Type Ia and II supernovae; reported interval is asymmetric.",
    },
    {
        "paper": "Cosmic flows in the nearby Universe: new peculiar velocities from SNe and cosmological constraints",
        "source_url": "https://doi.org/10.1093/mnras/staa2485",
        "velocity_km_s": 252,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.MULTI_TRACER,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 11,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "nearby universe",
        "effective_scale_mpc_h": "40",
        "year": 2020,
        "quality_score": 8,
        "notes": "Bulk flow at an effective radius of 40 h^-1 Mpc from a joint analysis including SNe Ia, SFI++, and 2MTF; not SNe-only.",
    },
    {
        "paper": "Peculiar-velocity cosmology with Types Ia and II supernovae (external flow beyond 2M++)",
        "source_url": "https://doi.org/10.1093/mnras/stab1446",
        "velocity_km_s": 195,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE_PREDICTED_FLOW_COMPARISON,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": (-23, 22),
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "low-z supernova sample",
        "effective_scale_mpc_h": "outside the 2M++ reconstruction volume",
        "year": 2021,
        "quality_score": 8,
        "notes": "Joint comparison of the 775-SN sample with the 2M++ reconstructed field; distinct from the 245 km/s SNe-only result in the same paper.",
    },
    {
        "paper": "A new measurement of the bulk flow of X-ray luminous clusters of galaxies (dark-flow claim)",
        "source_url": "https://doi.org/10.1088/2041-8205/712/1/L81",
        "velocity_km_s": None,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.CMB_KSZ,
        "result_kind": ResultKind.REPORTED_RANGE,
        "margin_of_error_km_s": None,
        "upper_limit_km_s": None,
        "reported_range_km_s": (600, 1000),
        "redshift_range": "claimed coherent to z ~ 0.3",
        "effective_scale_mpc_h": "~300-800 Mpc (reported physical scale)",
        "year": 2010,
        "quality_score": 2,
        "notes": "Controversial kSZ dark-flow claim; range is the authors' reported amplitude, not a Gaussian confidence interval. Later work has not established this as a consensus detection.",
    },
    {
        "paper": "Hawai'i Supernova Flows: Bulk Flow Measurements using SNe Ia in the Optical and NIR",
        "source_url": "https://arxiv.org/abs/2606.17142",
        "velocity_km_s": None,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.SNE,
        "result_kind": ResultKind.REPORTED_RANGE,
        "margin_of_error_km_s": None,
        "upper_limit_km_s": None,
        "reported_range_km_s": (100, 400),
        "redshift_range": "z ≲ 0.1",
        "effective_scale_mpc_h": None,
        "year": 2026,
        "quality_score": 7,
        "notes": "Recent preprint abstract summarizes results across samples/methods as roughly 100-400 km/s; it does not give one pooled estimate or uncertainty, so this is stored only as a reported range.",
    },
    {
        "paper": "Streaming motions of galaxy clusters within 12,000 km/s – V. The peculiar velocity field (Lauer & Postman result)",
        "source_url": "https://academic.oup.com/mnras/article/352/1/61/999349",
        "velocity_km_s": 689,
        "velocity_type": VelocityType.COSMOLOGICAL,
        "method": VelocityMethod.CLUSTER_DISTANCES,
        "result_kind": ResultKind.ESTIMATE,
        "margin_of_error_km_s": 178,
        "upper_limit_km_s": None,
        "reported_range_km_s": None,
        "redshift_range": "cluster sample to cz < 15,000 km/s",
        "effective_scale_mpc_h": "~80-110",
        "year": 1994,
        "quality_score": 5,
        "notes": "Historical brightest-cluster-galaxy distance result; later surveys and analyses have reported substantially different flows.",
    },
)


def velocityDataFrame():
    """Return the curated literature table of velocity measurements.

    ``velocity_km_s`` and all errors/limits are in km/s. A symmetric margin of
    error is a scalar; an asymmetric interval is ``(lower, upper)``; separate
    statistical/systematic errors are a mapping. Upper limits and author-
    reported ranges are stored in their own columns rather than misrepresented
    as point estimates. Distance scales are in h^-1 Mpc unless stated in notes.
    ``year`` is publication year. ``quality_score`` is a subjective 1-10
    appraisal, not a source-paper rating.
    """
    return pd.DataFrame.from_records(_MEASUREMENTS)


__all__ = ["VelocityType", "VelocityMethod", "ResultKind", "velocityDataFrame"]
