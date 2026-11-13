#=
Compare the update methods of FixedPointToolkit on a set of canonical fixed-point problems.

Run from the repository root:

    julia --project=. Sample/MixingComparison.jl            # table only
    julia --project=. Sample/MixingComparison.jl --plot     # also save convergence plots (needs Plots in your environment)

What to look at: the number of F calls. In real self-consistency loops (mean-field theory, DFT, DMFT, ...)
one call of F costs milliseconds to hours, so the method that needs the fewest calls wins, and the time
spent inside the update itself is irrelevant.

Rules of thumb that come out of this table:
- SimpleMixing: fine when F is cheap and the iteration is already a good contraction; crawls on slow modes.
- AndersonMixing: the default choice for expensive F. O(N m) memory; m = 5-10 is usually enough for
  nonlinear problems, while stiff linear problems want a longer history.
- BroydenMixing: very robust but keeps a dense N×N matrix, so use it when N is at most a few thousand.
- Newton's method is not included: a finite-difference Jacobian costs N + 1 calls of F per step
  (367 calls on the N = 60 problem below, versus about 35 for Anderson).
=#

using FixedPointToolkit, LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "test", "problems.jl"))

const TOL       =   1e-8
const MAX_ITER  =   2000
const PLOT      =   "--plot" in ARGS

if PLOT
    @eval using Plots
    plotdir = mkpath(joinpath(@__DIR__, "plots"))
end

##  Best of 3 runs after a compiling warm-up run.
function timed_run(p, Update, kwargs)
    sc, ncalls = run_problem(p, Update, kwargs; tol = TOL, max_iter = MAX_ITER)
    t = minimum(@elapsed(run_problem(p, Update, kwargs; tol = TOL, max_iter = MAX_ITER)) for _ in 1:3)
    return sc, ncalls, t
end

methods  = mixing_methods()
problems = canonical_problems()

println("F calls to reach max|F(x) - x| < $(TOL)   (- = not converged in $(MAX_ITER) iterations; time in ms)\n")
@printf("%-30s", "problem")
for (mname, _) in methods
    @printf("%17s", mname)
end
println()

for p in problems
    @printf("%-30s", p.name)
    PLOT && (plt = plot(yscale = :log10, xlabel = "F calls", ylabel = "max |F(x) - x|", title = p.name, titlefontsize = 10))

    for (mname, (Update, kwargs)) in methods
        sc, ncalls, t = timed_run(p, Update, kwargs)
        converged = residual(p, sc.VIns[end]) < TOL
        @printf("%17s", converged ? @sprintf("%d (%.1f)", ncalls, 1e3 * t) : @sprintf("- (%.1f)", 1e3 * t))

        if PLOT
            history = [maximum(abs, o .- i) for (i, o) in zip(sc.VIns, sc.VOuts)]
            plot!(plt, max.(history, 1e-16), label = mname, lw = 2)
        end
    end
    println()
    PLOT && savefig(plt, joinpath(plotdir, replace(p.name, r"[^A-Za-z0-9]+" => "_") * ".png"))
end

PLOT && println("\nPlots saved in $(plotdir)")
