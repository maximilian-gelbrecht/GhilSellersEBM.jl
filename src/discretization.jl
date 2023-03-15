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
struct Grid{T, V<:AbstractVector{T}, I<:Integer}
    N::I
    x::V
    Δx::T
end 

function Grid(x::AbstractRange{T}) where T
    N = length(x)
    Δx = abs(x[2] - x[1])
    return Grid(N, x, Δx)
end

"""
    NeumannFD{T}

1D - Finite Difference Scheme matrix with Neumann boundary conditions

# Initialization

    NeumannFD(T::DataType, n::Integer, Δx::Number=1, order="4th")

* `T`: precision used (e.g. `Float64`)
* `n`: number of grid points 
* `Δx`: spacing between grid points
* `order`: uses either a '4th'-order scheme or a '2nd'-order scheme, '4th' order scheme uses '2nd' order scheme at the points neighbouring the boundaries

    NeumannFD(grid::Grid{T})

* `grid`: instance of [`Grid`](@ref)

# Usage 

```julia
g = Grid(1:1:10)
∂ₓ = NeumannFD(g)
dfdx = ∂ₓ(f)
```
"""
struct NeumannFD{T,A<:AbstractArray{T,2}} 
    M::A
end 

function NeumannFD(T::DataType, n::Integer, Δx::Number=1, order="4th")  
    if order == "2nd"
        M = diagm(-1=>(-1*ones(T, n-1)),1=>ones(T, n-1))
        M[1,2] = T(0)
        M[n,n-1] = T(0)
        M ./= T(2*Δx)
        return NeumannFD(sparse(M))
    else 
        order == "4th"
        M = diagm(-1=>[T(-1)/T(2); (-(T(2)/T(3)).*ones(T, n-4)); T(-1)/T(2); T(0)],1=>[T(0); T(1)/T(2); (T(2)/T(3)).*ones(T, n-4); T(1)/T(2)])
        M += diagm(-2=>[((T(1)/T(12)).*ones(T, n-4)); T(0); T(0)],2=>[T(0); T(0); -(T(1)/T(12)).*ones(T, n-4)])
        M ./= T(Δx)
        return NeumannFD(sparse(M))
    end
end 

(FD::NeumannFD{T})(x::V) where {T,V<:AbstractVector{T}} = FD.M * x

NeumannFD(grid::Grid{T}, order="4th") where T = NeumannFD(T, grid.N, grid.Δx, order)