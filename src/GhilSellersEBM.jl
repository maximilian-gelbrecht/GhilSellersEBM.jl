module GhilSellersEBM

using LinearAlgebra, Interpolations, SparseArrays

include("discretization.jl")
include("parameters.jl")
include("model.jl")

export ghilsellers_ebm!, Grid, ContinousGhilSellersParameters

end
