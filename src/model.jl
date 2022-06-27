"""
    R_i(T,p::ContinousGhilSellersParameters)

Incoming radiation at discretized coordinate `x` and temperature `T` 
"""
R_i(T,p::ContinousGhilSellersParameters) = p.μ .* p.Q .* (1 .- α(T,p))

"""
    α(T,p::ContinousGhilSellersParameters)

Albedo at discretized coordinate `x` and temperature `T`. Albedo is cutoff at minimum 0.25 and maximum 0.85
"""
α(T,p::ContinousGhilSellersParameters) = cutoff_function.(p.b - p.c_1 .* (p.T_m .+ min.(0, T - p.c_2 .* p.z .- p.T_m)))


function cutoff_function(x) 
    if x <= 0.25
        return 0.25
    elseif x >= 0.65 
        return 0.65
    else 
        return x
    end 
end

"""
    R_o(T,p::ContinousGhilSellersParameters)

Outgoing radiation according to Stefan-Boltzmann law with emissivity reduced to account for greenhouse gases
"""
R_o(T,p::ContinousGhilSellersParameters) = p.σ .* T.^4 .* c(T,p)  

"""
    c(T,p::ContinousGhilSellersParameters)

Emissivity coefficient accounting for greenhouse gases
"""
c(T,p::ContinousGhilSellersParameters) = 1 .- p.m .* tanh.(p.c_3 .* T.^6)

"""
    k(T,p::ContinousGhilSellersParameters)

Heat flux, first term sensible heat flux, second term latent heat flux
"""
k(T,p::ContinousGhilSellersParameters) = p.k_1 + p.k_2 .* g(T,p)
k(T,i::Integer,p::ContinousGhilSellersParameters) = p.k_1[i] + p.k_2[i] * g(T,p)

"""
    g(T,p::ContinousGhilSellersParameters)

Latent heat flux contribution
"""
g(T,p::ContinousGhilSellersParameters) = p.c_4 .* exp.(-p.c_5 ./ T) ./ (T.^2)


"""
    D(T,p::ContinousGhilSellersParameters)  

1D Diffusion with Neumann boundary conditions applied
"""
function D(T,p::ContinousGhilSellersParameters)  
    A = ((2/π)^2) .* (1 ./ cos.(p.ϕ))
    
    D = A .* p.∂ₓ( cos.(p.ϕ) .* k(T,p) .* p.∂ₓ(T))
    D[1] = 2 * A[1] * cos(p.ϕ[1]) * k(T[1],1,p) * (T[2] - T[1])/(p.g.Δx^2)
    D[end] = 2 * A[end] * cos(p.ϕ[end]) * k(T[end],p.g.N,p) * (T[end-1] - T[end])/(p.g.Δx^2)

    return D
end

"""
    ghilsellers_ebm!(du,u,p,t)   

RHS of the PDE, to be used with `DifferentialEquations.jl`: `du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C`

# Usage 

```julia
using GhilSellersEBM, DifferentialEquations

x = (-90.:5.:90.)./90.
grid = Grid(x)
p = ContinousGhilSellersParameters(grid);
tspan = (0.,1e9)
prob = ODEProblem(ghilsellers_ebm!, 280*ones(p.g.N), tspan, p)
sol = solve(prob)
````
"""
function ghilsellers_ebm!(du,u,p,t)   
    du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C  
end