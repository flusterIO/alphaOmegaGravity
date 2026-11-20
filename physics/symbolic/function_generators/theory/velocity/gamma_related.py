from sympy import I, sqrt
from physics.symbolic.variables import beta


def gammaImaginary():
    r"""
    This is derived from the relationship in the `model_equation_reference` file.
    $$
    \arcsin{\left( x \right)} = -i \ln \left( i x + \sqrt{1 - x^{2}} \right)
    $$

    Then since $\gamma$ is the integral of arcsin, differentiating this equation for x gives $\gamma$.
    """
    return -I * (-beta / sqrt(1 - beta**2) + I) / (I * beta + sqrt(1 - beta**2))
