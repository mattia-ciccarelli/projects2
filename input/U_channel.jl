function load_input()

file="U_channel.msh";

# MATERIAL: density and dyn viscosity
Re = 5200 # bifurcation at 4140
Uin = 1
# coeffnorm = 100.0
coeffnorm = 100.0
material = [[1.0, coeffnorm/Re]]

# SOLID: associates materials to sets
solid = [[6, 1]]

# nodal DBC: each row a bc. node, direction, val
dbc = [[1,2],[2,1],[2,2],[4,2]]
dbcval = [0.,0.,0.,0.]

# inflow region
inflowbc = [[5,1],[5,2]]

# define parameter for the analysis
info=Sinfo()
info.neig=200   # number of modes to be computed
info.Lmm = [201, 288]  # master modes
info.Lmmcj = [567, 654]  # and their conjugates
info.Ffreq=1  # mode number that will give the freq (only one but +iomegat and -iomegat)
#info.Fmodes=[1]  # loading will be prop to sum of these modes
#info.Fmult=[0.5]  # with these amplitudes
info.style = 'c'
info.max_order = 3
info.max_orderNA = 0
info.tol=1e-1

#info.LambdaSym = [-1, -1, 1, 1, 0]
#info.LambdaSym = [-1, 1, 0]
info.LambdaSym = [-1, -1, 1, 1, 0]

return info, file, Re, Uin, coeffnorm, material, solid, dbc, dbcval, inflowbc

end
