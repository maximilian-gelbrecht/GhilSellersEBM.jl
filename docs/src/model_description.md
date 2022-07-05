# Model Description 

The model and its implementation is based on the following papers: 

* [Sellers: "A Global Climatic Model Based on the Energy Balance of the Earth-Atmosphere System", 1969](https://journals.ametsoc.org/view/journals/apme/8/3/1520-0450_1969_008_0392_agcmbo_2_0_co_2.xml)
* [Ghil: "Climate Stability for a Sellers-Type Model", 1976](https://journals.ametsoc.org/view/journals/atsc/33/1/1520-0469_1976_033_0003_csfast_2_0_co_2.xml)
* [Bodai et al: "Global instability in the Ghil-Sellers model", 2014](https://arxiv.org/abs/1402.3269)
    
What follows is a step by step introduction to the implementation. 

The model is a 1d energy balance model with the basic ansatz:

$$\begin{align} 
C(x) \frac{\partial T}{\partial t} &= R_{i}(x,T) - R_{o}(x,T) + D(x,T,\partial T/\partial x, \partial^2 T/\partial x^2)
\end{align}$$ 

with 
* the dimensionless rescaled latitude $x = 2\phi/\pi$
* the latitude dependent heat capacity $C(x)$ 
* the in- and outgoing radiation  $R_i$ and $R_o$ 
* the 1d diffusion  $D$ 
* Neumann boundary conditions $\frac{\partial T}{\partial x}(-1,t) = \frac{\partial T}{\partial x}(1,t) = 0\forall t$ are applied

## Discretization 

The model is discretized with a second-order central finite difference scheme. For this purpose [`Grid`](@ref) and [`NeumannFD`](@ref) are defined 

```@docs
GhilSellersEBM.Grid 
GhilSellersEBM.NeumannFD
```

## Parametriziation 

* The model includes a lot of physical parameters that were estimated from observational data to account for processes that are not resolved in the model itself. 
* For example it is not possible to directly model how the effective heat capacity $C(x)$ depends on latitude with a simple model, so the heat capacity at different latitudes is estimated from data and used a parameter in the model
* Sellers' paper gives these parameters as tabled values at two different grid, each with 10 degree spacing but offset from each other 
* We have to have all parameters on the same grid and we might want to solve the model with different grid spacings as well
* So we will take the values from the paper and do a cubic interpolation so that we can have parameter values at all grid points we need
* For this purpose we use `Interpolations.jl` and its `CubicSplineInterpolation`
* The parameters are symmetric with respect to the hemisphere
* Corrections suggested by [Bodai et al](https://arxiv.org/abs/1402.3269) are used

The Julia package provide [`ContinousGhilSellersParameters`](@ref) for this: 

```@docs
ContinousGhilSellersParameters
```

## Incoming Radiation 

* The incoming radiation is very similar to those of 0D-EBMs. It features an empirically estimated temperature and latitude dependent albedo 
* The solar irradiance is also latitude dependend and one of the parameters

$$\begin{align}
R_i(x,T) = \mu Q(x)(1 - \alpha(x,T))\end{align}$$

* The form of the albedo is emperically found: 

$$\begin{align}\alpha(x,T) = \left\{b(x) - c_1[T_m + \min(T - c_2 z(x) - T_m,0)]\right\}_c\\
&\text{with}\quad \{ x \}_c = \begin{cases} 0.25\quad\text{if } x<= 0.25\\ 
x\quad\text{if }0.25<x<0.65\\   0.65\quad\text{if } x>= 0.65 \end{cases}  \end{align}$$

* It is bound by the cutoff function $\{ x \}_c$ to be at least $0.25$ (so a completely icefree Earth) but at maximum $0.65$ (a completely frozen over Earth). 
* Ghil and Sellers used a maximum of $0.85$, Bodai et al reported that $0.65$ leads to results that are more in line with more comprehensive models 
* $\mu$ is an adjustable parameter to produce a modulate the solar irradiance in order to produce a bifurcation diagram. $\mu=1$ corresponds to present day conditions

```@docs 
GhilSellersEBM.R_i
GhilSellersEBM.α
```

## Outgoing Radiation 

* The outgoing radiation follows Stefan-Boltzmann law 
* A temperature dependend emissivity coefficient is introduced that reduces the outgoing radiation when temperature increases as greenhouse gases increase

$$\begin{align}
R_o(T) = \sigma c(T) T^4\end{align}$$

```@docs 
GhilSellersEBM.R_o
GhilSellersEBM.c
```

## Diffusion 

* The diffusion is given by

$$\begin{align} 
D(x,T) = \left(\frac{2}{\pi}\right)^2 \frac{1}{\cos(\phi(x))}\frac{\partial}{\partial x}\left(\cos(\phi(x))k(x,T) \frac{\partial T}{\partial x}\right)
\end{align}$$
#### Boundary Conditions

For the diffusion we have to make sure that the Neumann boundary conditions are fullfilled: 

$$ \frac{\partial T}{\partial x}(-1,t) = 0\quad \frac{\partial T}{\partial x}(1,t) = 0 \quad\forall t$$ 

For this purpose we have to rewrite the equation a bit, with $A(x)= ((2/π)^2) * (1 /\cos(ϕ(x)))$ and $B(x,T)=\cos(\phi(x)) k(x,T)$

$$\begin{align}
\frac{\partial T}{\partial t} &= R_i(x,T) - R_o(x,T) + A \frac{\partial}{\partial x}\left( B(x,T) \frac{\partial}{\partial x}T\right) \\
&= R_i(x,T) - R_o(x,T) + A \frac{\partial B}{\partial x}\frac{\partial T}{\partial x} + AB \frac{\partial^2 T}{\partial x^2} \\ 
\end{align}$$

Now, we only the only derivatives with respect to $x$ in the equation are $\frac{\partial T}{\partial x}$. So, we can look at the boundary condition: 

$$\begin{align} \frac{\partial T}{\partial x}(-1,t) &= 0 \\
\frac{T_2 - T_0}{2\Delta x} &= 0\quad\rightarrow\quad T_2 = T_0 \end{align}$$

Let's plug that into the equations above, to get the dynamics for the boundary points: 

$$\begin{align}\frac{\partial T_1}{\partial t} &= R_i(T_1) - R_o(T_1) + A \frac{\partial B}{\partial x} \left(\frac{T_2 - T_0}{2\Delta x}\right) + AB\frac{T_2 - 2T_1 + T_0}{\Delta x^2} \\
&= R_i(T_1) - R_o(T_1) + 2AB\frac{T_2 - T_1}{\Delta x^2}\\ 
\frac{\partial T_N}{\partial t} &= R_i(T_N) - R_o(T_N) + 2AB\frac{T_{N-1} - T_N}{\Delta x^2}\end{align}$$

```@docs
GhilSellersEBM.k
GhilSellersEBM.g
GhilSellersEBM.D
```
## Putting Everything Together 

* We now have everything for the model!
* Either use the predefined [`ghilsellers_ebm!`](@ref) or set up your own equation that you can modify

$$\begin{align} 
\frac{\partial T}{\partial t} &= (R_{i}(x,T) - R_{o}(x,T) + D(x,T,\partial T/\partial x, \partial^2 T/\partial x^2))/C(x)
\end{align}$$ 

```julia 
function ghilsellers_ebm!(du,u,p,t)   
    du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C  
end
```

```@docs
GhilSellersEBM.ghilsellers_ebm!    
```