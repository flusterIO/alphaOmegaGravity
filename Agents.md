# Agent guide

## IMPORTANT

Export fucking everything from `./physics/symbolic/__init__.py`, `./physics/plot_generators/__init__.py`, `./physics/symbolic/function_generators/electromagnetism/four_potential/__init__.py`, or `./physics/functions/__init__.py` according to whether the function is primarily symbolic, numeric, or plot related. Don't export shit from the nested `__init__.py` files..., and for fucks sake, dont' fucking remove any existing imports that are already there!

Use ONLY the tunnel files listed above. DO NOT EXPORT ANYTHING FROM ANYWHERE ELSE, AND DON"T MOVE ANY FUCKING IMPORTS.

## Project overview

This repository contains the `physics` Python package. It includes numerical and symbolic physics functions, constants and values, geometry helpers, tensor utilities, plotting helpers, and notebook support. Preserve the existing public API and the mathematical conventions used by the surrounding code when making changes.

## Repository layout

- `physics/symbolic/variables.py` defines shared SymPy coordinate systems, symbols, functions, and substitution data.
- `physics/symbolic/function_generators/` contains symbolic physics functions, grouped by subject area.
- `physics/symbolic/operators/`, `physics/symbolic/tensors_kinda/`, and `physics/symbolic/evaluators/` contain symbolic operators, tensor helpers, and evaluation utilities.
- `physics/functions/` contains the package's non-symbolic function implementations.
- `physics/values/`, `physics/data/`, and `physics/bodies/` contain value, data, and physical-body support code.
- `physics/plot_generators/` and `physics/notebook/` contain plotting and notebook helpers.

## Import rules

- **Always import shared symbolic variables directly from `physics.symbolic.variables` to avoid circular imports.** For example, use `from physics.symbolic.variables import Sx, Sy, Sz, c, phi`.
- Do not import shared variables through `physics.symbolic` from code inside the `physics.symbolic` package. Its `__init__.py` imports many function generators, so importing through the package can create circular dependencies.
- Import only the names needed by a module instead of relying on wildcard imports.
- Keep imports at module scope unless a local import is specifically needed to break a dependency cycle; do not use local imports to work around the direct-variable-import rule.

## Symbolic physics conventions

- Use SymPy expressions and matrices for symbolic results; do not coerce symbolic values to Python numbers or NumPy arrays unless the function's purpose requires numerical evaluation.
- Reuse the shared coordinates and symbols from `physics.symbolic.variables` so expressions compose and compare consistently across modules.
- Preserve the coordinate order and component order used by the function's stated mathematical convention. For four-vectors and matrices, document the convention in the function's docstring when it is not obvious.
- Use `smp` as the SymPy alias where that matches the surrounding modules.
- Keep functions focused and deterministic: construct and return the symbolic value described by the function, without unrelated substitutions or side effects.
- Be careful with units, sign conventions, index placement, and factors such as `c`; do not change them without checking the governing formula and nearby usage.

## Code changes

- Follow the existing Python style in the nearest files, including naming and formatting conventions. The repository contains both older camelCase names and newer snake_case names; preserve established names and signatures when editing existing APIs.
- Keep changes scoped to the requested behavior. Avoid unrelated cleanup or generated-file changes.
- Add or update docstrings when a function's purpose, return shape, units, or coordinate convention would otherwise be unclear.
- Do not edit generated artifacts such as `__pycache__` files.

## Validation

- When asked to validate a change, prefer a small import or focused check that exercises the changed behavior, followed by relevant project checks if their configuration is available.
- Report exactly which checks were run and whether they passed. Do not claim a test or check was run when it was not.
