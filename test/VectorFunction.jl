using LinearAlgebra, Test, FixedPointToolkit, Random
"""
Testing of the self consistency code for simpler multivariate functions using all implemented methods
"""

function f(v::Vector{Float64}) :: Vector{Float64}
    x = v[1]
    y = v[2]

    return [(x^2+y^2+8)/10, (x*y^2+x+8)/10]
end

Random.seed!(1234)
guess                  =   rand(2)
mixingAlpha            =   0.5 

selfcons               =   SelfCons(f, SimpleMixing   , guess; Update_kwargs = Dict(:alpha => mixingAlpha))
selfcons_scheduled     =   SelfCons(f, ScheduledMixing, guess; Update_kwargs = Dict(:alpha => mixingAlpha, :Scheduler => ExponentialScheduling))
selfcons_broyden       =   SelfCons(f, BroydenMixing  , guess; Update_kwargs = Dict(:alpha => mixingAlpha, :B => mixingAlpha * Matrix(1.0I, length(guess), length(guess)) ))
selfcons_anderson      =   SelfCons(f, AndersonMixing , guess; Update_kwargs = Dict(:alpha => mixingAlpha, :m => 5))

tol                    =   1e-8
FixedPoint!(selfcons          ; tol = tol, max_iter=100)
FixedPoint!(selfcons_scheduled; tol = tol, max_iter=100)
FixedPoint!(selfcons_broyden  ; tol = tol, max_iter=100)
FixedPoint!(selfcons_anderson ; tol = tol, max_iter=100)

@test norm(f(selfcons.VIns[end]) - selfcons.VIns[end]) < tol * sqrt(length(guess))
@test norm(f(selfcons_scheduled.VIns[end]) - selfcons_scheduled.VIns[end]) < tol * sqrt(length(guess))
@test norm(f(selfcons_broyden.VIns[end]) - selfcons_broyden.VIns[end]) < tol * sqrt(length(guess))
@test norm(f(selfcons_anderson.VIns[end]) - selfcons_anderson.VIns[end]) < tol * sqrt(length(guess))
