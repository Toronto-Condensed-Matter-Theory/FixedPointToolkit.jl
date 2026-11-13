# FixedPointToolkit.jl

FixedPointToolkit.jl is a Julia package meant for calculating fixed points of any scalar or vector functions, or equivalently, solve some system of self-consistent equations.

Currently supported :
* Fixed point iteration can work on any user defined function.
* User has options to use different pre-defined update methods to be used during the iterations : Simple Mixing, Scheduled Mixing, Broyden Mixing, and Anderson (Pulay/DIIS) Mixing. See `Sample/MixingComparison.jl` for a comparison on standard test problems.
* Can also implement a custom update function for user-specific scenarios.
* Can allow for different kinds of scheduling where the step-size in parameter space during the iterations changes. Right now supports Quadratic and Exponential scheduling, but can also support custom user-defined scheduler as well.
* Can checkpoint and save results into JLD2 files, and resume iterations from reading such files.
* Can plot results of inputs, outputs, and convergence as a function of iterations (load `Plots` first: `using Plots`).

