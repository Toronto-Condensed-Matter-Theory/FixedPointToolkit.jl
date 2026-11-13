##  Interrupt a run, resume it from the JLD2 checkpoint, and check the mixing state survives exactly:
##  the resumed run must take the same number of F calls as an uninterrupted one.

for (mname, (Update, kwargs)) in mixing_methods()
    occursin(r"Broyden|Anderson", mname) || continue
    p = problem_mixed_spectrum()

    @testset "$(mname)" begin
        mktempdir() do dir
            file = joinpath(dir, "checkpoint.jld2")

            full, nfull = run_problem(p, Update, kwargs; tol = 1e-8, max_iter = 500)

            sc = with_logger(NullLogger()) do
                sc = SelfCons(p.F, Update, p.x0; Update_kwargs = copy(kwargs))
                FixedPoint!(sc, file; tol = 1e-8, max_iter = 5, checkpoint_interval = 1)
                sc
            end

            saved = read_checkpoint(file)
            state = occursin("Broyden", mname) ? (:B,) : (:X, :G)
            for key in state
                @test saved["Update kwargs"][key] == sc.Update_kwargs[key]
            end

            resumed = with_logger(NullLogger()) do
                ContinueFixedPoint!(file, p.F, Update, Dict{Symbol, Real}(:max_iter => 500, :tol => 1e-8, :checkpoint_interval => 100))
            end
            @test residual(p, resumed.VIns[end]) < 1e-8
            @test length(resumed.VOuts) == nfull
        end
    end
end
