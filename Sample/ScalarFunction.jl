using FixedPointToolkit

using LinearAlgebra

"""
Lets find the fixed point of f(x) = cos(x) using three different methods.
"""

guess                  =   rand()
mixingAlpha            =   0.5   

selfcons               =   SelfCons(cos, SimpleMixing   , guess; Update_kwargs = Dict(:alpha => mixingAlpha))
selfcons_scheduled     =   SelfCons(cos, ScheduledMixing, guess; Update_kwargs = Dict(:alpha => mixingAlpha, :Scheduler => QuadraticScheduling))
selfcons_broyden       =   SelfCons(cos, BroydenMixing  , guess; Update_kwargs = Dict(:alpha => mixingAlpha, :B => mixingAlpha * Matrix(1.0I, length(guess), length(guess)) ))

FixedPoint!(selfcons          ; tol=1e-8, max_iter=100)
FixedPoint!(selfcons_scheduled; tol=1e-8, max_iter=100)
FixedPoint!(selfcons_broyden  ; tol=1e-8, max_iter=100)

