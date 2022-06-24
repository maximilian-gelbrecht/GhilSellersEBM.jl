using GhilSellersEBM
using Test, OrdinaryDiffEq, StatsBase

@testset "GhilSellersEBM.jl" begin
    # Write your tests here.

    x = (-90.:5.:90.)./90.
    grid = Grid(x)
    p = ContinousGhilSellersParameters(grid);

    include("fd.jl")
    include("integration_test.jl")
end

