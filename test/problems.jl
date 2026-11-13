#=
Canonical fixed-point test problems, shared by the test suite and `Sample/MixingComparison.jl`.

Every problem is a NamedTuple
    (name, F, x0, F_args, F_kwargs, exact)
where `F(x, F_args...; F_kwargs...)` is the map whose fixed point x = F(x) is sought, `x0` the
starting guess and `exact` the known fixed point (or `nothing`).
=#

using LinearAlgebra, Random, Logging

"""P1: x = cos(x). The fixed point is the Dottie number."""
function problem_cos()
    dottie = 0.7390851332151607
    return (name = "cos(x) (scalar)", F = cos, x0 = 0.3, F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = dottie)
end

"""P2: a small polynomial system with fixed point (1, 1)."""
function problem_poly2d()
    F(v) = [(v[1]^2 + v[2]^2 + 8) / 10, (v[1] * v[2]^2 + v[1] + 8) / 10]
    return (name = "2D polynomial", F = F, x0 = [0.2, 0.7], F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = [1.0, 1.0])
end

"""
P3: Chandrasekhar H-equation, discretised with the composite midpoint rule on N nodes,

    H(μ_i) = 1 / (1 - (c/2) Σ_j  μ_i H_j / (μ_i + μ_j) / N).

The classic benchmark for Anderson acceleration (Kelley); it becomes hard as c → 1.
"""
function problem_chandrasekhar(; N::Int = 100, c::Float64 = 0.99)
    μ   = ((1:N) .- 0.5) ./ N
    K   = (c / (2N)) .* (μ ./ (μ .+ μ'))
    F(h; K = K) = 1 ./ (1 .- K * h)
    return (name = "Chandrasekhar H, c = $(c)", F = F, x0 = ones(N), F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = nothing)
end

"""
P4: Jacobi iteration for the 1D Poisson problem -u'' = 1 on N interior points.
Linear, with spectral radius cos(π/(N+1)): the slow, smooth modes make plain mixing crawl.
"""
function problem_jacobi_poisson(; N::Int = 30)
    h2  = 1 / (N + 1)^2
    function F(u)
        v = similar(u)
        for i in eachindex(u)
            left  = i > 1 ? u[i-1] : 0.0
            right = i < N ? u[i+1] : 0.0
            v[i]  = (h2 + left + right) / 2
        end
        return v
    end
    A     = SymTridiagonal(fill(2.0, N), fill(-1.0, N - 1))
    exact = A \ fill(h2, N)
    return (name = "Jacobi 1D Poisson, N = $(N)", F = F, x0 = zeros(N), F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = exact)
end

"""
P5: a nonlinear map x ↦ A tanh(x) + b whose linearisation has a spectrum spread over [-0.9, 0.97]:
both slow modes (λ → 1) and oscillating ones (λ → -1), like a frustrated self-consistency loop.
"""
function problem_mixed_spectrum(; N::Int = 60, seed::Int = 1)
    rng = Xoshiro(seed)
    Q   = Matrix(qr(randn(rng, N, N)).Q)
    A   = Q * Diagonal(collect(range(-0.9, 0.97, length = N))) * Q'
    b   = 0.3 .* randn(rng, N)
    x0  = 0.1 .* randn(rng, N)
    F(x) = A * tanh.(x) .+ b
    return (name = "Mixed spectrum tanh, N = $(N)", F = F, x0 = x0, F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = nothing)
end

"""
P6: a complex linear map x ↦ M x + c with complex eigenvalues of modulus up to 0.95.
Tests that the methods are correct for complex order parameters.
"""
function problem_complex_linear(; N::Int = 40, seed::Int = 2)
    rng = Xoshiro(seed)
    Q   = Matrix(qr(randn(rng, ComplexF64, N, N)).Q)
    λ   = range(0.3, 0.95, length = N) .* cis.(2π .* rand(rng, N))
    M   = Q * Diagonal(λ) * Q'
    c   = randn(rng, ComplexF64, N)
    F(x) = M * x .+ c
    return (name = "Complex linear, N = $(N)", F = F, x0 = zeros(ComplexF64, N), F_args = (), F_kwargs = Dict{Symbol, Any}(), exact = (I - M) \ c)
end

canonical_problems() = [
    problem_cos(),
    problem_poly2d(),
    problem_chandrasekhar(c = 0.9),
    problem_chandrasekhar(c = 0.99),
    problem_chandrasekhar(c = 0.9999),
    problem_jacobi_poisson(),
    problem_mixed_spectrum(),
    problem_complex_linear(),
]

"""Update methods to compare: name => (Update function, Update_kwargs)."""
mixing_methods() = [
    "Simple α=0.5"     => (SimpleMixing,    Dict{Symbol, Any}(:alpha => 0.5)),
    "Scheduled α=0.5"  => (ScheduledMixing, Dict{Symbol, Any}(:alpha => 0.5, :Scheduler => QuadraticScheduling)),
    "Broyden α=0.5"    => (BroydenMixing,   Dict{Symbol, Any}(:alpha => 0.5)),
    "Anderson m=5"     => (AndersonMixing,  Dict{Symbol, Any}(:alpha => 0.5, :m => 5)),
    "Anderson m=10"    => (AndersonMixing,  Dict{Symbol, Any}(:alpha => 0.5, :m => 10)),
]

"""Run one method on one problem; return the SelfCons and the number of F evaluations."""
function run_problem(p, Update, Update_kwargs; tol = 1e-8, max_iter = 2000)
    sc = with_logger(NullLogger()) do
        sc = SelfCons(p.F, Update, p.x0; F_args = p.F_args, F_kwargs = copy(p.F_kwargs), Update_kwargs = copy(Update_kwargs))
        FixedPoint!(sc; tol = tol, max_iter = max_iter)
        sc
    end
    return sc, length(sc.VOuts)     # one F call per stored output (incl. the initial one)
end

residual(p, x) = maximum(abs, p.F(x, p.F_args...; p.F_kwargs...) .- x)
