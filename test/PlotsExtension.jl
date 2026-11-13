ENV["GKSwstype"] = "100"     # headless GR
using Plots

@test Base.get_extension(FixedPointToolkit, :FixedPointToolkitPlotsExt) !== nothing

sc, _ = run_problem(problem_poly2d(), AndersonMixing, Dict{Symbol, Any}(:alpha => 0.5))
for arg in ("inputs", "outputs", "convergence")
    @test Plot_History(sc, arg) isa Plots.Plot
end
