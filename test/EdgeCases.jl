@testset "F_kwargs reach F in every update" begin
    ##  F *requires* the keyword `shift`, so a method that drops or mis-passes F_kwargs errors out.
    Fk(x; shift) = cos.(x) .+ shift
    for (mname, (Update, kwargs)) in mixing_methods()
        sc = with_logger(NullLogger()) do
            sc = SelfCons(Fk, Update, [0.1, 0.5]; F_kwargs = Dict{Symbol, Any}(:shift => 0.1), Update_kwargs = copy(kwargs))
            FixedPoint!(sc; tol = 1e-10, max_iter = 500)
            sc
        end
        @test maximum(abs, Fk(sc.VIns[end]; shift = 0.1) .- sc.VIns[end]) < 1e-10
    end
end

@testset "Broyden at the fixed point keeps a finite B" begin
    p   = problem_poly2d()
    x   = p.exact
    out = BroydenMixing(x, p.F(x), p.F)
    @test all(isfinite, out["kwargs"][:B])
    @test out["VInNext"] ≈ x
end

@testset "Broyden secant condition" begin
    p   = problem_mixed_spectrum()
    x   = p.x0
    out = BroydenMixing(x, p.F(x), p.F; alpha = 0.5)
    B   = out["kwargs"][:B]
    dx  = out["VInNext"] - x
    dg  = (out["VOutNext"] - out["VInNext"]) - (p.F(x) - x)
    @test B * dg ≈ -dx
end

@testset "Anderson history is bounded" begin
    p  = problem_mixed_spectrum()
    for m in (0, 1, 3)
        kw = Dict{Symbol, Any}(:alpha => 0.5, :m => m)
        x, y = p.x0, p.F(p.x0)
        for _ in 1:10
            out  = AndersonMixing(x, y, p.F; kw...)
            kw   = out["kwargs"]
            x, y = out["VInNext"], out["VOutNext"]
            @test length(kw[:X]) == length(kw[:G]) <= m + 1
        end
    end
    ##  Scalars: at most one difference can ever be independent.
    out = AndersonMixing(0.3, cos(0.3), cos; m = 5, X = [[0.1], [0.2]], G = [[cos(0.1) - 0.1], [cos(0.2) - 0.2]])
    @test length(out["kwargs"][:X]) == 2
end

@testset "Anderson with m = 0 is simple mixing" begin
    p  = problem_mixed_spectrum()
    sA, nA = run_problem(p, AndersonMixing, Dict{Symbol, Any}(:alpha => 0.5, :m => 0))
    sS, nS = run_problem(p, SimpleMixing,   Dict{Symbol, Any}(:alpha => 0.5))
    @test nA == nS
    @test all(isapprox.(sA.VIns, sS.VIns; rtol = 1e-10))
end

@testset "Anderson drops an ill-conditioned history" begin
    ##  Two identical history entries make ΔG singular; the step must fall back to simple mixing, not NaN.
    x   = [0.2, 0.4]
    g   = cos.(x) .- x
    out = AndersonMixing(x, cos.(x), x -> cos.(x); alpha = 0.5, m = 3, X = [x], G = [g])
    @test out["VInNext"] ≈ x .+ 0.5 .* g
    @test length(out["kwargs"][:X]) == 1
end

@testset "Anderson restarts when the residual grows" begin
    x   =   [0.2, 0.4]
    F   =   v -> cos.(v)
    g   =   F(x) .- x
    ##### a previous residual much smaller than the current one: history must be cleared (simple-mixing step)
    out =   AndersonMixing(x, F(x), F; alpha = 0.5, m = 3, X = [x .+ 0.1], G = [1e-3 .* g])
    @test out["VInNext"] ≈ x .+ 0.5 .* g
    @test length(out["kwargs"][:X]) == 1
    ##### restart = Inf keeps the history
    out =   AndersonMixing(x, F(x), F; alpha = 0.5, m = 3, X = [x .+ 0.1], G = [1e-3 .* g], restart = Inf)
    @test length(out["kwargs"][:X]) == 2
end
