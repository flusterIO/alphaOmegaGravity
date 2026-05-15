# Included Functions

| Derivation                                                | Function                         | Error                              |
| --------------------------------------------------------- | -------------------------------- | ---------------------------------- |
| cosmological velocity from local gravitational parameters | `vbarCosmological`               | Difficult to assess                |
| kinematic velocity from local gravitational parameters    | `vbarKinematic`                  | $0.43\% - 0.92\%$                  |
| Earth's rotational velocity via $\alpha$                  | `rotationVelocityFromDivergence` | $0.006778\%$                       |
| $\alpha$, $2-\alpha$, $1 \pm \alpha$                      | `alpha`                          | $0.025\%$ <br/> $0.0034\%$ (solar) |
| Earth's cosmological velocity as a function of $\alpha$   | `vbarAlphaEquivalence`           | $0.08\%$                           |

> Don't forget the `describe` function that generates a similar note at runtime, using CODATA and IAU values from the Scipy and Astropy packages.

## The lighthouse & the Clocktower

> An opening thought experiment

- `vbarCosmological`

This arises from $\gamma$ and the equivalence principle. If $g$ is the result of spatial dilation instead of time dilation, then this equation is the velocity required to produce $g$ on Earth's surface.

$$
\bar{v}_0^{\prime} = c \sqrt{1 - \frac{1}{\left( 1 + \frac{g}{R} \right)^{2}}} = 525.49 ~\text{km}~\text{s}
$$

> This aligns closely with observation correlating distance SNe data to the CMB, but this is difficult to assess due to the large error margins inherent in these types of measurements.

Notice here that the units are in $km~s$, aka: distance.

Once time has been accounted for by assigning some equivalent quantities to distance[^1], we can derive a time dependent equivalent function as:

- `vbarKinematic`

$$
\bar{v}_0 = c \sqrt{1 - \frac{1}{\left( 1 + \frac{1}{2}\frac{g}{R} x \right)^{2}}} = 371.56 ~\text{km}~\text{s}^{-1}
$$

> This accurate relative to our CMB to well within 1% of observation.

The above equations in turn produce a dilation of space on Earth's surface equivalent to $g$ per unit _'time'_.

Since $\vec{g}$ is in it's linear form, this gives a vector, so taking it's derivative gives us some component representing the component of a 3-vec. This gives some function where $d\vec{s} = \frac{d}{dR} \vec{g}$.

This inverse of this determinant then gives $\omega$.

$$
\omega = \sqrt[3]{2\frac{GM}{R^3}} = 0.145377 \frac{m^{1/3}}{s^{1/3}}
$$

where

$$
\nabla \cdot \mathbf{S} = \left|\begin{array}{ccc}
     \omega & 0 & 0 \\
    0 &  \omega & 0 \\
    0 & 0 &  \omega
    \end{array}\right| = \omega^3 = \delta
$$

This is almost precisely $2 \alpha$, but we can do better. If we inspect the spinor function approximation, we find:

$$
\alpha =  \Gamma \sqrt[3]{\frac{GM}{R^3}}
$$

where

- `spinorApproximation`

$$
\Gamma \approx \frac{1}{2 - \alpha}
$$

This gives a combined equation of:

$$
\alpha = \frac{1}{2 - \alpha} \sqrt[3]{\frac{d}{dR}g} \quad \pm 0.025\%
$$

- `alphaSolar`

## Electromagnetism

- `rotationalVelocityFromDivergence`

Now let's account for the fact that if we're in co-motion along our temporal axis, then this spatial dilation would remain hidden.

Our experiments would measure $\nabla \cdot B = 0$, but in reality we're measuring $ \nabla \cdot B = \delta$. It only remains hidden due to the dilation of our measuring apparatus _with_ the dilation of space. If this is true, then

$$
\nabla \cdot E = c \left( \nabla \cdot B \right) = \delta c = 921.11628 \frac{m^2}{s^2}
$$

If we apply our spinor approximation, we find:

$$
\left(\frac{1}{2- \alpha}\right) \nabla \cdot \vec{E} = \frac{1}{2- \alpha} \delta c = 462.24472 \frac{m^2}{s^2} = \frac{2 \pi R_{\oplus}}{T_{\text{day}}} \pm 0.006778\%
$$

Implying that $\alpha$ must be in units of $m/s$, because that's Earth's rotational velocity to within $0.006778\%$. Another 'pseudo-scientific', 'dimensional coincidence' I guess?

> Note here that our polar radius was used to calculate this value to this degree of accuracy, as to eliminate any influence from centripital force. The result for our equatorial radius gives $0.34\%$ error.

## Our Velocity

- `vbarAlphaEquivalence`

In something like a geometric transformation between the dimensionalities of linear motion and divergence, $\alpha$ can be found as a core function producing our cosmological velocity. If $\alpha$ can be quantized, so too then can _motion_ and gravity.

$$
 \beta^2 = \frac{v_{\text{cosmological}}^2}{c^2} = \left( \alpha \left( 2 - \alpha \right) \right)^3 \pm 0.08\%
$$

## Mass

Consider the classical mass density formula after substituting $ \delta$ we find:

$$
\rho = \frac{3 \left( \frac{\delta R^{3}}{2G} \right)}{4 \pi R^{3}} = \frac{3 \delta}{8 \pi G}
$$

Since we found that the dimensionless ratio $ds/dx = \delta v_{0} \approx \Phi$, then

$$
\rho = \frac{\Phi}{\left( \frac{8\pi}{3}G \cdot v_{0} \right)}
$$

Using the constraint we found above, that $\beta = \left( \alpha \left( 2 - \alpha \right) \right)^{2/3}$

$$
\rho = \frac{\delta}{G}\cdot \frac{\Phi/\left( c \delta \right)}{\alpha \left( 2 - \alpha \right)}
$$

Or

$$
\rho = \frac{\Phi}{Gc \cdot \alpha \left( 2 - \alpha \right)}
$$

[^1]: See the section on 'time equivalent quantities' for a full derivation, but it's comically rudimentary.
