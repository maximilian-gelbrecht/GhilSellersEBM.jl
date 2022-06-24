"""
    Grid{T}

1D discretization Grid 

# Initialization 

Grid(x::AbstractRange{T})

* `x`: coordinate axis 

# Fields

* `N`: Number of grid points 
* `x`: Coordinate axis
* `Δx`: Spacing between grid points 
"""
struct Grid{T}
N 
x::AbstractVector{T}
Δx::T
end 

function Grid(x::AbstractRange{T}) where T
N = length(x)
Δx = abs(x[2] - x[1])
return Grid{T}(N, x, Δx)
end

"""
    NeumannFD{T}

1D - Finite Difference Scheme matrix with Neumann boundary conditions

# Initialization

    NeumannFD(T::DataType, n::Integer, Δx::Number=1)

* `T`: precision used (e.g. `Float64`)
* `n`: number of grid points 
* `Δx`: spacing between grid points

    NeumannFD(grid::Grid{T})

# Usage 

```julia
g = Grid(1:1:10)
∂ₓ = NeumannFD(g)
dfdx = ∂ₓ(f)
```
"""
struct NeumannFD{T} 
    M::AbstractArray{T,2}
end 

function NeumannFD(T::DataType, n::Integer, Δx::Number=1)
    M = diagm(-1=>(-1*ones(T, n-1)),1=>ones(T, n-1))
    M[1,2] = T(0)
    M[n,n-1] = T(0)
    M ./= T(2*Δx)
    NeumannFD(M)
end 

(FD::NeumannFD{T})(x::AbstractVector{T}) where T = FD.M * x

NeumannFD(grid::Grid{T}) where T = NeumannFD(T, grid.N, grid.Δx)