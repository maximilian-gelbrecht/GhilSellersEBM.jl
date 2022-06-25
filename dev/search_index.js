var documenterSearchIndex = {"docs":
[{"location":"model_description/#Model-Description","page":"Model Description","title":"Model Description","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The model and its implementation is based on the following papers: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Sellers: \"A Global Climatic Model Based on the Energy Balance of the Earth-Atmosphere System\", 1969\nGhil: \"Climate Stability for a Sellers-Type Model\", 1976\nBodai et al: \"Global instability in the Ghil-Sellers model\", 2014","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"What follows is a step by step introduction to the implementation. ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The model is a 1d energy balance model with the basic ansatz:","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign \nC(x) fracpartial Tpartial t = R_i(xT) - R_o(xT) + D(xTpartial Tpartial x partial^2 Tpartial x^2)\nendalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"with ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"the dimensionless rescaled latitude x = 2phipi\nthe latitude dependent heat capacity C(x) \nthe in- and outgoing radiation  R_i and R_o \nthe 1d diffusion  D \nNeumann boundary conditions fracpartial Tpartial x(-1t) = fracpartial Tpartial x(1t) = 0forall t are applied","category":"page"},{"location":"model_description/#Discretization","page":"Model Description","title":"Discretization","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The model is discretized with a second-order central finite difference scheme. For this purpose Grid and NeumannFD are defined ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"GhilSellersEBM.Grid \nGhilSellersEBM.NeumannFD","category":"page"},{"location":"model_description/#GhilSellersEBM.Grid","page":"Model Description","title":"GhilSellersEBM.Grid","text":"Grid{T}\n\n1D discretization Grid \n\nInitialization\n\nGrid(x::AbstractRange{T})\n\nx: coordinate axis \n\nFields\n\nN: Number of grid points \nx: Coordinate axis\nΔx: Spacing between grid points \n\n\n\n\n\n","category":"type"},{"location":"model_description/#GhilSellersEBM.NeumannFD","page":"Model Description","title":"GhilSellersEBM.NeumannFD","text":"NeumannFD{T}\n\n1D - Finite Difference Scheme matrix with Neumann boundary conditions\n\nInitialization\n\nNeumannFD(T::DataType, n::Integer, Δx::Number=1)\n\nT: precision used (e.g. Float64)\nn: number of grid points \nΔx: spacing between grid points\nNeumannFD(grid::Grid{T})\ngrid: instance of Grid\n\nUsage\n\ng = Grid(1:1:10)\n∂ₓ = NeumannFD(g)\ndfdx = ∂ₓ(f)\n\n\n\n\n\n","category":"type"},{"location":"model_description/#Parametriziation","page":"Model Description","title":"Parametriziation","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The model includes a lot of physical parameters that were estimated from observational data to account for processes that are not resolved in the model itself. \nFor example it is not possible to directly model how the effective heat capacity C(x) depends on latitude with a simple model, so the heat capacity at different latitudes is estimated from data and used a parameter in the model\nSellers' paper gives these parameters as tabled values at two different grid, each with 10 degree spacing but offset from each other \nWe have to have all parameters on the same grid and we might want to solve the model with different grid spacings as well\nSo we will take the values from the paper and do a cubic interpolation so that we can have parameter values at all grid points we need\nFor this purpose we use Interpolations.jl and its CubicSplineInterpolation\nThe parameters are symmetric with respect to the hemisphere\nCorrections suggested by Bodai et al are used","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The Julia package provide ContinousGhilSellersParameters for this: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"ContinousGhilSellersParameters","category":"page"},{"location":"model_description/#GhilSellersEBM.ContinousGhilSellersParameters","page":"Model Description","title":"GhilSellersEBM.ContinousGhilSellersParameters","text":"ContinousGhilSellersParameters{T}\n\nHolds all parameters needed for the Ghil Sellers 1D EBM. The model is given in CGS units.  Uses the tabled values from the Ghil Paper and interpolates on a new grid. \n\nInitialization\n\nContinousGhilSellersParameters(g::Grid, μ=1.)\n\ng::Grid{T} \nϕ::AbstractVector{T} latitude vector\n∂ₓ::NeumannFD{T} FD scheme\nT_0::AbstractVector{T} initial condition given in the Ghil paper\nμ::T modifier for the incoming solar irradiance, can be used to produce a bifuraction diagram\nC::AbstractVector{T} heat capacity\nQ::AbstractVector{T} solar irradiance\nb::AbstractVector{T} empirical albedo coefficient\nz::AbstractVector{T} average surface elevation \nk_1::AbstractVector{T} sensible heat flux coefficient\nk_2::AbstractVector{T}coeffcient of latent heat flux in the atmosphere \nc_1::T = 0.009 empirical albedo coefficient\nc_2::T = 0.0065 empirical albedo coefficient\nc_3::T = 1.9e-15 empirical emissivity coeffcient\nc_4::T = 6.105*0.75*exp(19.6) latent heat flux coefficient \nc_5::T = 5350 latent heat flux coefficient \nσ::T = 1.356e-12 Stefan Boltzmann constant (in CGS system)\nm::T = 0.5 Atmospheric attenuation coefficient (0.5 for present conditions)\nT_m::T = 283.16 reference temperature \n\n\n\n\n\n","category":"type"},{"location":"model_description/#Incoming-Radiation","page":"Model Description","title":"Incoming Radiation","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The incoming radiation is very similar to those of 0D-EBMs. It features an empirically estimated temperature and latitude dependent albedo \nThe solar irradiance is also latitude dependend and one of the parameters","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign\nR_i(xT) = Q(x)(1 - alpha(xT))endalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The form of the albedo is emperically found: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalignalpha(xT) = leftb(x) - c_1T_m + min(T - c_2 z(x) - T_m0)right_c\ntextwithquad  x _c = begincases 025quadtextif  x= 025 \nxquadtextif 025x065   065quadtextif  x= 065 endcases  endalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"It is bound by the cutoff function  x _c to be at least 025 (so a completely icefree Earth) but at maximum 065 (a completely frozen over Earth). \nGhil and Sellers used a maximum of 085, Bodai et al reported that 065 leads to results that are more in line with more comprehensive models ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"GhilSellersEBM.R_i\nGhilSellersEBM.α","category":"page"},{"location":"model_description/#GhilSellersEBM.R_i","page":"Model Description","title":"GhilSellersEBM.R_i","text":"R_i(T,p::ContinousGhilSellersParameters)\n\nIncoming radiation at discretized coordinate x and temperature T \n\n\n\n\n\n","category":"function"},{"location":"model_description/#GhilSellersEBM.α","page":"Model Description","title":"GhilSellersEBM.α","text":"α(T,p::ContinousGhilSellersParameters)\n\nAlbedo at discretized coordinate x and temperature T. Albedo is cutoff at minimum 0.25 and maximum 0.85\n\n\n\n\n\n","category":"function"},{"location":"model_description/#Outgoing-Radiation","page":"Model Description","title":"Outgoing Radiation","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The outgoing radiation follows Stefan-Boltzmann law \nA temperature dependend emissivity coefficient is introduced that reduces the outgoing radiation when temperature increases as greenhouse gases increase","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign\nR_o(T) = sigma c(T) T^4endalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"GhilSellersEBM.R_o\nGhilSellersEBM.c","category":"page"},{"location":"model_description/#GhilSellersEBM.R_o","page":"Model Description","title":"GhilSellersEBM.R_o","text":"R_o(T,p::ContinousGhilSellersParameters)\n\nOutgoing radiation according to Stefan-Boltzmann law with emissivity reduced to account for greenhouse gases\n\n\n\n\n\n","category":"function"},{"location":"model_description/#GhilSellersEBM.c","page":"Model Description","title":"GhilSellersEBM.c","text":"c(T,p::ContinousGhilSellersParameters)\n\nEmissivity coefficient accounting for greenhouse gases\n\n\n\n\n\n","category":"function"},{"location":"model_description/#Diffusion","page":"Model Description","title":"Diffusion","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"The diffusion is given by","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign \nD(xT) = left(frac2piright)^2 frac1cos(phi(x))fracpartialpartial xleft(cos(phi(x))k(xT) fracpartial Tpartial xright)\nendalign","category":"page"},{"location":"model_description/#Boundary-Conditions","page":"Model Description","title":"Boundary Conditions","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"For the diffusion we have to make sure that the Neumann boundary conditions are fullfilled: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"$","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"\\frac{\\partial T}{\\partial x}(-1,t) = 0\\quad \\frac{\\partial T}{\\partial x}(1,t) = 0 \\quad\\forall t$ ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"For this purpose we have to rewrite the equation a bit, with A(x)= ((2π)^2) * (1 cos(ϕ(x))) and B(xT)=cos(phi(x)) k(xT)","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign\nfracpartial Tpartial t = R_i(xT) - R_o(xT) + A fracpartialpartial xleft( B(xT) fracpartialpartial xTright) \n= R_i(xT) - R_o(xT) + A fracpartial Bpartial xfracpartial Tpartial x + AB fracpartial^2 Tpartial x^2  \n = R_i(xT) - R_o(xT) + A fracpartial Bpartial Tfracpartial Tpartial xfracpartial Tpartial x + AB fracpartial^2 Tpartial x^2  \n = R_i(xT) - R_o(xT) + A fracpartial Bpartial Tleft(fracpartial Tpartial xright)^2+ AB fracpartial^2 Tpartial x^2 \nendalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Now, we only the only derivatives with respect to x in the equation are fracpartial Tpartial x. So, we can look at the boundary condition: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign fracpartial Tpartial x(-1t) = 0 \nfracT_2 - T_02Delta x = 0quadrightarrowquad T_2 = T_0 endalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"Let's plug that into the equations above, to get the dynamics for the boundary points: ","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalignfracpartial T_1partial t = R_i(T_1) - R_o(T_1) + A fracpartial Bpartial T left(fracT_2 - T_02Delta xright)^2 + ABfracT_2 - 2T_1 + T_0Delta x^2 \n= R_i(T_1) - R_o(T_1) + 2ABfracT_2 - T_1Delta x^2 \nfracpartial T_Npartial t = R_i(T_N) - R_o(T_N) + 2ABfracT_N-1 - T_NDelta x^2endalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"GhilSellersEBM.k\nGhilSellersEBM.g\nGhilSellersEBM.D","category":"page"},{"location":"model_description/#GhilSellersEBM.k","page":"Model Description","title":"GhilSellersEBM.k","text":"k(T,p::ContinousGhilSellersParameters)\n\nHeat flux, first term sensible heat flux, second term latent heat flux\n\n\n\n\n\n","category":"function"},{"location":"model_description/#GhilSellersEBM.g","page":"Model Description","title":"GhilSellersEBM.g","text":"g(T,p::ContinousGhilSellersParameters)\n\nLatent heat flux contribution\n\n\n\n\n\n","category":"function"},{"location":"model_description/#GhilSellersEBM.D","page":"Model Description","title":"GhilSellersEBM.D","text":"D(T,p::ContinousGhilSellersParameters)\n\n1D Diffusion with Neumann boundary conditions applied\n\n\n\n\n\n","category":"function"},{"location":"model_description/#Putting-Everything-Together","page":"Model Description","title":"Putting Everything Together","text":"","category":"section"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"We now have everything for the model!\nEither use the predefined ghilsellers_ebm! or set up your own equation that you can modify","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"beginalign \nfracpartial Tpartial t = (R_i(xT) - R_o(xT) + D(xTpartial Tpartial x partial^2 Tpartial x^2))C(x)\nendalign","category":"page"},{"location":"model_description/","page":"Model Description","title":"Model Description","text":"function ghilsellers_ebm!(du,u,p,t)   \n    du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C  \nend","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GhilSellersEBM.Grid\nGhilSellersEBM.NeumannFD\nGhilSellersEBM.ContinousGhilSellersParameters\nGhilSellersEBM.R_i\nGhilSellersEBM.α\nGhilSellersEBM.R_o\nGhilSellersEBM.c\nGhilSellersEBM.k\nGhilSellersEBM.g\nGhilSellersEBM.D\nGhilSellersEBM.ghilsellers_ebm!    ","category":"page"},{"location":"reference/#GhilSellersEBM.ghilsellers_ebm!","page":"Reference","title":"GhilSellersEBM.ghilsellers_ebm!","text":"ghilsellers_ebm!(du,u,p,t)\n\nRHS of the PDE, to be used with DifferentialEquations.jl: du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C\n\nUsage\n\n```julia using GhilSellersEBM, DifferentialEquations\n\nx = (-90.:5.:90.)./90. grid = Grid(x) p = ContinousGhilSellersParameters(grid); tspan = (0.,1e9) prob = ODEProblem(ghilsellers_ebm!, 280*ones(p.g.N), tspan, p) sol = solve(prob) ````\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GhilSellersEBM","category":"page"},{"location":"#GhilSellersEBM","page":"Home","title":"GhilSellersEBM","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GhilSellersEBM. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package provides an implementation of the Ghil Sellers 1D Energy Balance Model. It is based on three publications: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Sellers: \"A Global Climatic Model Based on the Energy Balance of the Earth-Atmosphere System\", 1969\nGhil: \"Climate Stability for a Sellers-Type Model\", 1976\nBodai et al: \"Global instability in the Ghil-Sellers model\", 2014","category":"page"},{"location":"#Example-Use","page":"Home","title":"Example Use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Solving the model with two initial conditions: one leading to cold and one to a warm state. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GhilSellersEBM, DifferentialEquations, Plots \n\nx = (-90.:5.:90.)./90.\ngrid = Grid(x)\np = ContinousGhilSellersParameters(grid);\n\ntspan = (0.,1e9)\nprob = ODEProblem(ghilsellers_ebm!, 220*ones(p.g.N), tspan, p)\n\nsol_1 = solve(prob)\nsol_2 = solve(remake(prob, u0=290*ones(p.g.N)))\n\nt_plot = range(tspan[1],tspan[2],length=200)\nanim = @animate for it ∈ t_plot\n    plot(p.ϕ, sol_1(it), xlabel=\"Latitude ϕ [rad]\", label=\"Cold State\", ylims=[210,300], ylabel=\"Temperature T [K]\", title=\"Temperature Profile of Ghil Sellers EBM\")\n    plot!(p.ϕ, sol_2(it), label=\"Warm State\", ylims=[200,300])\nend \ngif(anim, \"ebm-anim.gif\", fps=10)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: ../figures/ebm-anim.gif)","category":"page"}]
}
