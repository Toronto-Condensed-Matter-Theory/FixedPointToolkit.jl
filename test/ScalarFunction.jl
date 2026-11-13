using Test, FixedPointToolkit, LinearAlgebra, Random

"""
Lets find the fixed point of f(x) = cos(x) using the different methods.
"""

Random.seed!(1234)
guess                  =   rand()
mixingAlpha            =   0.5   

selfcons               =   SelfCons(cos, SimpleMixing   , guess; Update_kwargs = Dict(:alpha => mixingAlpha))
selfcons_scheduled     =   SelfCons(cos, ScheduledMixing, guess; Update_kwargs = Dict(:alpha => mixingAlpha, :Scheduler => QuadraticScheduling))
selfcons_broyden       =   SelfCons(cos, BroydenMixing  , guess; Update_kwargs = Dict(:alpha => mixingAlpha, :B => mixingAlpha * Matrix(1.0I, 1, 1) ))
selfcons_anderson      =   SelfCons(cos, AndersonMixing , guess; Update_kwargs = Dict(:alpha => mixingAlpha, :m => 3))

tol                    =   1e-8

FixedPoint!(selfcons          ; tol = tol, max_iter = 100)
FixedPoint!(selfcons_scheduled; tol = tol, max_iter = 100)
FixedPoint!(selfcons_broyden  ; tol = tol, max_iter = 100)
FixedPoint!(selfcons_anderson ; tol = tol, max_iter = 100)

@test norm(cos(selfcons.VIns[end]) - selfcons.VIns[end]) < tol
@test norm(cos(selfcons_scheduled.VIns[end]) - selfcons_scheduled.VIns[end]) < tol
@test norm(cos(selfcons_broyden.VIns[end]) - selfcons_broyden.VIns[end]) < tol
@test norm(cos(selfcons_anderson.VIns[end]) - selfcons_anderson.VIns[end]) < tol
@test selfcons_anderson.VIns[end] isa Float64
