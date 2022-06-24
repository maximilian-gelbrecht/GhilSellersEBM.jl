# test that the finite difference scheme works 
g2 = Grid(range(0,2π,step=π/10)) # initialize the grid 
cos_data = cos.(g2.x) # evaluate cos on the grid 
∂_x = NeumannFD(g2) # initialize the FD 
dcos_data = ∂_x(cos_data) # apply the FD 
sin_data = -1 .* sin.(g2.x)

@test sum(abs2,sin_data -  dcos_data) < 0.01
