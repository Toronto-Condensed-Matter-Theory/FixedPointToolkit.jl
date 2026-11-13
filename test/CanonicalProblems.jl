##  Every method on every canonical problem (see problems.jl), plus iteration-count expectations.

const TOL       =   1e-8
const MAX_ITER  =   2000

##  Problems where plain (and scheduled) mixing cannot reach TOL within MAX_ITER.
too_slow_for_mixing(name) = occursin("Poisson", name) || occursin("0.9999", name)

##  Measured F-call counts (Julia 1.13) for (Broyden, Anderson m=5), used as regression caps with 2× headroom.
const MEASURED  =   Dict(
    "cos(x) (scalar)"               => (6, 6),
    "2D polynomial"                 => (10, 8),
    "Chandrasekhar H, c = 0.9"      => (12, 13),
    "Chandrasekhar H, c = 0.99"     => (18, 14),
    "Chandrasekhar H, c = 0.9999"   => (34, 23),
    "Jacobi 1D Poisson, N = 30"     => (36, 232),
    "Mixed spectrum tanh, N = 60"   => (45, 38),
    "Complex linear, N = 40"        => (45, 71),
)

calls = Dict{Tuple{String, String}, Int}()

for p in canonical_problems(), (mname, (Update, kwargs)) in mixing_methods()
    occursin(r"Simple|Scheduled", mname) && too_slow_for_mixing(p.name) && continue

    sc, ncalls          =   run_problem(p, Update, kwargs; tol = TOL, max_iter = MAX_ITER)
    calls[(p.name, mname)] = ncalls
    x                   =   sc.VIns[end]

    @testset "$(mname) on $(p.name)" begin
        @test residual(p, x) < TOL
        @test eltype(x) == eltype(p.x0)
        if p.exact !== nothing
            @test maximum(abs, x .- p.exact) < 1e-5
        end
    end
end

@testset "Accelerated methods beat simple mixing" begin
    for p in canonical_problems()
        bcap, acap = 2 .* MEASURED[p.name]
        @test calls[(p.name, "Broyden α=0.5")] <= bcap
        @test calls[(p.name, "Anderson m=5")]  <= acap

        if !too_slow_for_mixing(p.name) && !occursin(r"cos|polynomial|0\.9,", p.name)
            for m in ("Broyden α=0.5", "Anderson m=5", "Anderson m=10")
                @test calls[(p.name, m)] < calls[(p.name, "Simple α=0.5")]
            end
        end
    end

    ##  Regression for the pre-1.2 sign convention, which needed 250 calls here vs 132 for simple mixing.
    @test calls[("Mixed spectrum tanh, N = 60", "Broyden α=0.5")] < calls[("Mixed spectrum tanh, N = 60", "Simple α=0.5")] ÷ 2
end
