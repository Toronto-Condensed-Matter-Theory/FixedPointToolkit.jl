using FixedPointToolkit
using Test, LinearAlgebra, Random, Logging

include("problems.jl")

@testset "FixedPointToolkit.jl" begin
    @testset "Scalar function"      include("ScalarFunction.jl")
    @testset "Vector function"      include("VectorFunction.jl")
    @testset "Canonical problems"   include("CanonicalProblems.jl")
    @testset "Update edge cases"    include("EdgeCases.jl")
    @testset "Checkpointing"        include("Checkpointing.jl")
    @testset "Plots extension"      include("PlotsExtension.jl")
end
