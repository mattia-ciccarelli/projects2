
file="cylinder.msh";

# MATERIAL: density and dyn viscosity
Re = 51.0
# Re=49.03266832152 #bifurcation
# Re = 20.0
Uin = 1.5
coeffnorm = 100.0
# coeffnorm = 1.0
material = [[1.0, coeffnorm*0.1/Re]]

# SOLID: associates materials to sets
solid = [[3, 1]]

# nodal DBC: each row a bc. node, direction, val
dbc = [[1,1],[1,2]]
dbcval = [0.,0.]

# inflow region
inflowbc = [[2,1],[2,2]]

# define parameter for the analysis
info=Sinfo()
info.neig=2    # number of modes to be computed
info.Lmm = [1]  # master modes
info.Lmmcj = [2]  # and their conjugates
info.Ffreq=1  # mode number that will give the freq (only one but +iomegat and -iomegat)
#info.Fmodes=[1]  # loading will be prop to sum of these modes
#info.Fmult=[0.5]  # with these amplitudes
info.style = 'c'
info.max_order = 3
info.max_orderNA = 0
info.tol=1e-1

# info.LambdaSym = [-1, -1, 1, 1, 0]
info.LambdaSym = [-1, 1, 0]


