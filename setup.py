from setuptools import setup, find_packages

setup(
    name="alpha_omega_gravity",
    packages=find_packages(),
    version="0.1.4",
    description="""Successfully unified gravity and electromagnetism.""",
    long_description_content_type="text/markdown",
    long_description="""
# $\\alpha \\omega$ Gravity

> Classically simultaneous General Relativity via an unwrapped Kaluza-Klein dimension

Yes, that's $\\alpha$, as in the fine-structure constant. $\\omega$ is defined within the package, because I can't resist the urge to be grandiose. 


## What is this?

This is my personal toolkit for the development of the theory of gravity described [here](https://zenodo.org/records/24280888). It might not work flawlessly as I've hardcoded some path's for my own use case as default parameters, but if the paper takes off and you want to see how it came to be, feel free to explore the code base. The perfect workflow for me was just a Ghostty/PtPython repl session, the plugin available [here](https://github.com/flusterIO/symbolic_math_repl), and this package.

## What's the story behind this?

In the very early weeks of 2021 I was putting away my work for the night at around 3am when I put in some toast, grabbed a book I had sitting next to me, and flipped to a random page. The book was a collection of Einstein's lectures and papers, and the random line I flipped to was his famous quote about the length of the moving rod from "On The Electrodynamics Of Moving Bodies".

Maybe it was just that it was out of context, I've read that line a dozen times before, but it occurred to me... Einstein's assumptions made far more sense before the observations that gave us the Big Bang model. I spent the next few weeks playing with the math before I realized that absolutely everything aligns with observation when you bind the Lorentz factor to gravitational acceleration and cosmic inflation, so I quit my job. I slept on the side of the road for 3 years and 4 months to the day, but I think I have something significant.

If you found this package for some reason before you find the paper, you can take a read [here](https://zenodo.org/records/24280888).

    """,
)
