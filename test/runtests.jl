using GhilSellersEBM
using Test, OrdinaryDiffEq, StatsBase

@testset "GhilSellersEBM.jl" begin
    include("fd.jl")
    include("integration_test.jl")
end

