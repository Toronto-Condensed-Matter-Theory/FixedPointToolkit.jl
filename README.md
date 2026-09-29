# FixedPointToolkit.jl

[![Build Status](https://github.com/Toronto-Condensed-Matter-Theory/FixedPointToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Toronto-Condensed-Matter-Theory/FixedPointToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

FixedPointToolkit.jl is a Julia package meant for calculating fixed points of any scalar or vector functions, or equivalently, solve some system of self-consistent equations.

Currently supported :
* Fixed point iteration can work on any user defined function.
* User has options to use different pre-defined update methods to be used during the iterations : Simple Mixing, Schduled Mixing, and Broyden Mixing.
* Can also implement a custom update function for user-specific scenarios.
* Can allow for different kinds of scheduling where the step-size in parameter space during the iterations changes. Right now supports Quadratic and Exponential scheduling, but can also support custom user-defined scheduler as well.
* Can checkpoint and save results into JLD2 files, and resume iterations from reading such files.
* Can plot results of inputs, outputs, and convergence as a function of iterations.

## Documentation
For further details, please refer to [Documentation](https://Toronto-Condensed-Matter-Theory.github.io/FixedPointToolkit.jl/dev/)
