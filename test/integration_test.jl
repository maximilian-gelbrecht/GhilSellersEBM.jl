# test that the EBM integrates to a cold and warm state
x = (-90.:5.:90.)./90.
grid = Grid(x)
p = ContinousGhilSellersParameters(grid);

tspan = (0.,1e9)
prob = ODEProblem(ghilsellers_ebm!, 220*ones(p.g.N), tspan, p)

# cold state 
sol_1 = solve(prob, Tsit5())
@test 215 <= mean(sol_1) <= 225

# warm state
sol_2 = solve(remake(prob, u0=280*ones(p.g.N)), Tsit5())
@test 270 <= mean(sol_2) <= 285