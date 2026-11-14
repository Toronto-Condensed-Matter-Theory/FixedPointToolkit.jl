#=
Regenerate the README figures of FixedPointToolkit.jl.

    julia --project=docs/readme docs/readme/figures.jl

Every figure comes from FixedPointToolkit runs on the canonical test problems in `test/problems.jl`.
The figures are written to `docs/src/assets/readme/`.
=#

include("style.jl")
using FixedPointToolkit, LinearAlgebra, Logging
include(joinpath(@__DIR__, "..", "..", "test", "problems.jl"))

const TOL       =   1e-8
const MAX_ITER  =   2000

##  One colour per update method, the same in every figure.
##  ScheduledMixing is left out: with a long iteration budget its step size stays close to alpha,
##  and its curve lies on top of SimpleMixing.
const METHODS = [
    (name = "Simple",   Update = SimpleMixing,   kwargs = Dict{Symbol, Any}(:alpha => 0.5),           color = COLORS.blue),
    (name = "Broyden",  Update = BroydenMixing,  kwargs = Dict{Symbol, Any}(:alpha => 0.5),           color = COLORS.orange),
    (name = "Anderson", Update = AndersonMixing, kwargs = Dict{Symbol, Any}(:alpha => 0.5, :m => 5),  color = COLORS.aqua),
]
const CALL_TICKS = ([10, 30, 100, 300, 1000], ["10", "30", "100", "300", "1000"])

##  Residual max|F(x) - x| after each call of F.
residual_history(sc) = [maximum(abs, o .- i) for (i, o) in zip(sc.VIns, sc.VOuts)]

##  Number of F calls needed to converge, or NaN if the run did not converge.
function calls_to_converge(p, Update, kwargs)
    sc, ncalls = run_problem(p, Update, kwargs; tol = TOL, max_iter = MAX_ITER)
    return residual(p, sc.VIns[end]) < TOL ? ncalls : NaN
end


##########  1. Convergence of the update methods on the Chandrasekhar H-equation  ##########
let p = problem_chandrasekhar(c = 0.99)
    plt = plot(yscale = :log10, xlabel = L"\mathrm{calls\ of\ } F", ylabel = L"\max_i\, |F(x)_i - x_i|",
               ylims = (1e-12, 10), yticks = 10.0 .^ (-12:2:0), legend = :topright)
    for m in METHODS
        sc, ncalls = run_problem(p, m.Update, m.kwargs; tol = TOL, max_iter = MAX_ITER)
        println(rpad(m.name, 10), ncalls, " calls")
        plot!(plt, max.(residual_history(sc), 1e-16), label = m.name, color = m.color)
    end
    hline!(plt, [TOL], color = COLORS.muted, linestyle = :dash, lw = 1, label = "")
    savefig_readme(plt, "convergence")
end


##########  2. Cobweb diagram of x = cos(x)  ##########
let
    sc  =   with_logger(NullLogger()) do
        sc = SelfCons(cos, SimpleMixing, 0.05; Update_kwargs = Dict{Symbol, Any}(:alpha => 1.0))
        FixedPoint!(sc; tol = 1e-6, max_iter = 100)
        sc
    end
    xs  =   range(0, 1.05, length = 400)
    plt =   plot(xs, cos.(xs), color = COLORS.ink, label = L"F(x) = \cos x",
                 xlabel = L"x", ylabel = L"F(x)", xlims = (0, 1.05), ylims = (0, 1.05),
                 aspect_ratio = :equal, size = (440, 420), legend = :bottomright)
    plot!(plt, xs, xs, color = COLORS.muted, linestyle = :dash, lw = 1.25, label = L"F(x) = x")

    ##  Staircase: (x_n, x_n) -> (x_n, F(x_n)) -> (F(x_n), F(x_n)) with x_{n+1} = F(x_n).
    cx  =   Float64[sc.VIns[1]];  cy = Float64[0.0]
    for (i, o) in zip(sc.VIns, sc.VOuts)
        push!(cx, i, o);  push!(cy, o, o)
    end
    plot!(plt, cx, cy, color = COLORS.blue, lw = 1.25, label = "iterates")
    scatter!(plt, [sc.VIns[end]], [sc.VOuts[end]], color = COLORS.orange, markersize = 6,
             markerstrokecolor = :white, label = "fixed point")
    println("cobweb: x* = ", sc.VIns[end], " after ", length(sc.VIns), " calls")
    savefig_readme(plt, "cobweb")
end


##########  3. Sensitivity to the mixing parameter alpha  ##########
let p = problem_mixed_spectrum()
    alphas  =   collect(0.05:0.05:1.0)
    plt     =   plot(yscale = :log10, xlabel = L"\mathrm{mixing\ parameter\ } \alpha", ylabel = L"\mathrm{calls\ of\ } F \mathrm{\ to\ converge}",
                     xlims = (0, 1.02), yticks = CALL_TICKS, legend = :topright)
    for m in METHODS
        calls = [calls_to_converge(p, m.Update, merge(m.kwargs, Dict{Symbol, Any}(:alpha => a))) for a in alphas]
        println(rpad(m.name, 10), Int.(replace(calls, NaN => -1)))
        plot!(plt, alphas, calls, label = m.name, color = m.color, marker = :circle, markersize = 4, markerstrokecolor = :white)
    end
    savefig_readme(plt, "alpha_scan")
end


##########  4. Cost as the Chandrasekhar problem approaches its singular point c -> 1  ##########
let
    gaps    =   10.0 .^ range(-1, -4, length = 13)
    plt     =   plot(xscale = :log10, yscale = :log10, xflip = true, xlims = (7.5e-5, 1.3e-1),
                     xticks = ([1e-1, 1e-2, 1e-3, 1e-4], [L"10^{-1}", L"10^{-2}", L"10^{-3}", L"10^{-4}"]), yticks = CALL_TICKS,
                     xlabel = L"1 - c", ylabel = L"\mathrm{calls\ of\ } F \mathrm{\ to\ converge}", legend = :topleft)
    for m in METHODS
        calls = [calls_to_converge(problem_chandrasekhar(c = 1 - g), m.Update, m.kwargs) for g in gaps]
        println(rpad(m.name, 10), Int.(replace(calls, NaN => -1)))
        plot!(plt, gaps, calls, label = m.name, color = m.color, marker = :circle, markersize = 4, markerstrokecolor = :white)
    end
    savefig_readme(plt, "chandrasekhar_scan")
end
