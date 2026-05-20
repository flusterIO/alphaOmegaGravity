# Divergent Gravity

### Results

| derivation                                                | function                         | % error             |
| :-------------------------------------------------------- | :------------------------------- | :------------------ |
| cosmological velocity from local gravitational parameters | `vbarCosmological`               | Difficult to assess |
| kinematic velocity from local gravitational parameters    | `vbarKinematic`                  | 0.442258391884371   |
| Earth's roational velocity via α                          | `rotationVelocityFromDivergence` | 0.00677756400547747 |
| Solar-Galactic Radius from visible mass                   | `galacticRadiusPolynomial`       | 0.0                 |
| α (solar)                                                 | `alphaSolar`                     | 0.00342749503676447 |
| α (Earth)                                                 | `omega`, `alphaFactor`           | 0.0253043303850287  |
| α (orbital velocity)                                      | `solarMetricVelocity`, `vOrbit`  | 0.0544665494810030  |
| α (from derived velocity)                                 | `vbarAlphaEquivalence`           | 0.0761820393770567  |
| J ≈ Φ                                                     | `currentEquivalent`              | 0.213633171310522   |
| h from α                                                  | `omega`, `alpha`                 | 0.0253887188634106  |
| h/(2 - α) from ω                                          | `omega`, `alpha`                 | 0.00342749503676447 |
| Φ mass density approximation                              | `massDensityEquivalence`         | 0.863214423082925   |

Well, if you're ever going to read a `readme.md` you picked the right one to read.

To make a long story short, I have a background in astrophysics but was working as a software engineer when I had a realization... Einstein's assumptions made far more sense before the observations that gave us the Big Bang.

I started to play around with the math, realized I could do something with it, became homeless, kept working on it, and 5 years later I think I solved it. $\alpha$, the fine-structure constant defined in terms of gravitational parameters to undeniable accuracy, _precisely_ where the model predicts to find it.

The problem? The model concludes that gravity is repulsive and that the Earth expands _with_ the rest of the Universe. Instead of celebrating, I'm called a crackpot, a 'nutjob', or any variety of names, but the math is valid. _This many_ alignments are not a 'dimensional coincidence', it's a real breakthrough in the way we understand gravity, whether or not people are willing to accept it.

This package will define $\alpha$ in terms of gravitational and velocity based parameters using a variety of bodies and functions, demonstrating that this is not just a one-off coincidence.

To go along with this model I built [Fluster](https://flusterapp.com), a pair of note taking applications for STEM students, with the intention of promoting this model. There is a more in-depth summary available with the application, and I'm hoping to have the online derivation up within a week (I'm still homeless, so finding wifi time is a challenge), but the entire model flows directly from a single, mathematically equivalent modification to Einstein's SR which then produces a direct physical mechanism to the equivalence principle. From there deriving our velocity relative to distant standard candles is rudimentary, but the complete derivation for our kinematic velocity, along with all of the derivations for the functions included here will be available online as soon as I get my websit back in working order.

But lasty, a closing thought experiment to justify the case for classical simultaneity:

### The Lighthouse & The Clocktower

> It is from this principle that all of the math in this package can be derived.

Consider two observers, $A$ and $B$. Let $A$ and $B$ agree that $B$ should travel between two arbitrary points in space[^2], $p_1$ and $p_2$ at an **_agreed_** upon velocity $v$.

At $p_2$ place a time keeping device that remains visible from $p_1$. Allow that $A$ remains at rest relative to $p_1$, and very near $p_1$ while $B$ goes in motion (at an inertial pace) between $p_1$ and $p_2$.

At what time does $A$ observe $B$ reach $p_2$? If the value is not $d/v$, but $\gamma d/v$, then $v \neq v_\prime$. The only possible solution that allows for Einstein's additional factor of $\gamma$ while allowing $A$ to remain in an unchanged state is that $d \to d\gamma$ in the coordinate system of the observer in relative motion, while $v \to v/\gamma$ in the coordinate system of the observer at relative rest to accomodate the change in spatial density.
