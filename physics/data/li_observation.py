import astropy.constants.iau2015 as ac
import numpy as np
import pandas as pd
import scipy.constants as spc

HEIGHTS = np.array(
    [
        5.62,
        7.0,
        9.87,
    ]
)

TIME_DELTA_AT_5_62 = 1.9 * 10**-12
TIME_DELTA_AT_7_0 = 2.7 * 10**-12
TIME_DELTA_AT_9_87 = 4.5 * 10**-12

TIME_DELAYS = np.array([TIME_DELTA_AT_5_62, TIME_DELTA_AT_7_0, TIME_DELTA_AT_9_87])

FIBER_LENGTH = 20_000

FIBER_GROUP_REFRACTIVE_INDEX = 1.4682
GROUP_VELOCITY_KM_S = 3e5 / FIBER_GROUP_REFRACTIVE_INDEX
GROUP_VELOCITY = GROUP_VELOCITY_KM_S * 1e3

ALPHA_LI_REPORTED = 1.37e-2

# Li's paper gives a comparison value derived from an
# Einstein gravitational-potential model.
ALPHA_EINSTEIN_COMPARISON = 1.6e-10  # s^2/m
# Instrument and acquisition details reported in the paper.
DELAY_MEASUREMENT_RESOLUTION = 0.1e-12  # s = 0.1 ps
MEASUREMENTS_PER_HEIGHT = 10  # repeated tests
DATA_POINTS_PER_HEIGHT = 1000  # more than 1000 reported
DATA_COLLECTION_DURATION = 30.0  # s per height
MEASUREMENT_INTERVAL = 120.0  # s between height comparisons

MEASURED_GRADIENT, FIT_INTERCEPT = np.polyfit(
    HEIGHTS,
    TIME_DELAYS,
    deg=1,
)

# Gradient in picoseconds per metre.
MEASURED_GRADIENT_PS_PER_M = MEASURED_GRADIENT * 1e12


FITTED_DELAYS = MEASURED_GRADIENT * HEIGHTS + FIT_INTERCEPT

RESIDUAL_SUM_OF_SQUARES = np.sum((TIME_DELAYS - FITTED_DELAYS) ** 2)

TOTAL_SUM_OF_SQUARES = np.sum((TIME_DELAYS - np.mean(TIME_DELAYS)) ** 2)

R_SQUARED = (
    1.0 - RESIDUAL_SUM_OF_SQUARES / TOTAL_SUM_OF_SQUARES
    if TOTAL_SUM_OF_SQUARES > 0
    else np.nan
)


ALPHA_FROM_GRAPH = (
    MEASURED_GRADIENT
    * ac.R_earth.value
    * GROUP_VELOCITY
    / (2.0 * spc.G * ac.M_earth.value / ac.R_earth.value**2 * FIBER_LENGTH)
)


PREDICTED_GRADIENT = (
    2.0
    * ALPHA_LI_REPORTED
    * spc.G
    * ac.M_earth.value
    / ac.R_earth.value**2
    * FIBER_LENGTH
    / (ac.R_earth.value * GROUP_VELOCITY)
)


PREDICTED_GRADIENT_PS_PER_M = PREDICTED_GRADIENT * 1e12

# Predicted delays at the three experimental heights.
PREDICTED_TIME_DELAYS = PREDICTED_GRADIENT * HEIGHTS

# Ratio of reported alpha to alpha inferred from the approximate
# points. A large discrepancy signals that the plotted values,
# parameter interpretation, or equation inputs need further review.
ALPHA_DISCREPANCY_RATIO = (
    ALPHA_LI_REPORTED / ALPHA_FROM_GRAPH if ALPHA_FROM_GRAPH != 0 else np.nan
)

data = np.array(
    [
        [5.62, TIME_DELTA_AT_5_62],
        [7.0, TIME_DELTA_AT_7_0],
        [9.87, TIME_DELTA_AT_9_87],
    ]
)


def liObservationDataFrame():
    """Return the per-height observational values and related delays.

    Heights are in metres, velocities in metres per second, and delays in
    seconds. The measured delay is the approximate value read from the paper's
    figure; fitted and predicted delays are included for comparison.
    """
    return pd.DataFrame(
        {
            "height_m": HEIGHTS,
            "group_velocity_m_s": np.full(HEIGHTS.shape, GROUP_VELOCITY),
            "measured_delay_s": TIME_DELAYS,
            "fitted_delay_s": FITTED_DELAYS,
            "paper_delta_t_s": (
                2.0
                * ALPHA_LI_REPORTED
                * (spc.G * ac.M_earth.value / ac.R_earth.value**2)
                * FIBER_LENGTH
                * HEIGHTS
                / (ac.R_earth.value * GROUP_VELOCITY)
            ),
        }
    )

if __name__ == "__main__":
    print("Li (2025) optical-fibre experiment")
    print("=" * 48)

    print(f"Fibre length: {FIBER_LENGTH:.1f} km")
    print(f"Group refractive index: {FIBER_GROUP_REFRACTIVE_INDEX}")
    print(f"Group velocity: {GROUP_VELOCITY:.6e} m/s")
    print(f"Delay measurement resolution: {DELAY_MEASUREMENT_RESOLUTION * 1e12:.1f} ps")

    print("\nApproximate data from Figure 2(b)")
    print("Height (m) | Delay (ps)")
    for height, delay in data:
        print(f"{height:10.2f} | {delay:9.3f}")

    print("\nDelay-versus-height fit")
    print(f"Gradient: {MEASURED_GRADIENT_PS_PER_M:.6f} ps/m")
    print(f"Gradient: {MEASURED_GRADIENT:.6e} s/m")
    print(f"Intercept: {FIT_INTERCEPT * 1e12:.6f} ps")
    print(f"R^2: {R_SQUARED:.6f}")

    print("\nModel parameter comparison")
    print(f"Li's reported alpha: {ALPHA_LI_REPORTED:.6e} s^2/m")
    print(f"Alpha from approximate graph: {ALPHA_FROM_GRAPH:.6e} s^2/m")
    print(f"Reported/inferred alpha ratio: {ALPHA_DISCREPANCY_RATIO:.3f}")

    print("\nGradient predicted by Li's reported alpha")
    print(f"Predicted gradient: {PREDICTED_GRADIENT_PS_PER_M:.6f} ps/m")

    print("\nPredicted delays using the reported alpha")
    for height, delay in zip(HEIGHTS, PREDICTED_TIME_DELAYS):
        print(f"{height:5.2f} m: {delay * 1e12:10.4f} ps")
