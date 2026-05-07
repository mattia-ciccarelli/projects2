
file="two_cylinders.msh";

# MATERIAL: density and kinematic viscosity
Re = 56                                     # Reynolds
rho = 1.0                                   # Density
d = 1.0                                     # Diameter of the obstacle
coeffnorm = 1.0                             # Normalisation for conditioning (haven't tested)
material = [[rho, coeffnorm*rho*d/Re]]
g_value = 0.74

# SOLID: associates materials to sets
solid = [[7, 1]]

# MESH MATERIAL: Young's modulus, Poisson's coefficient and density
# mesh_material = [[1, 0.499, 0.1]]
# mesh = [[7, 1]]

Uin = 1.0

# nodal DBC: each row a bc. node, direction, val
dbc = [[1,1],[1,2],[3,1],[3,2],[5,1],[5,2],[6,1],[6,2]]
dbcval = [Uin,0.0,Uin,0.0,0.0,0........................................0,0.0,0.0]

# INFLOW region
# ATTENTION: for this example the inflow velocity is constant, not parabolic!
#            Check if it's like this on the code.
inflowbc = [[4,1],[4,2]]

# fluid-structure interface
# fs_interface = [5,6]
# mesh_boundary = [1,2,3,4]

# define parameter for the analysis
info=Sinfo()
info.neig=20    # number of modes to be computed
info.Lmm = [1]#[28]#[31]#, 37]  # master modes
info.Lmmcj = [2]#[29]#[32]#, 38]  # and their conjugates
info.Ffreq=1  # mode number that will give the freq (only one but +iomegat and -iomegat)
#info.Fmodes=[1]  # loading will be prop to sum of these modes
#info.Fmult=[0.5]  # with these amplitudes
info.style = 'c'
info.max_order = 3
info.max_orderNA = 0
info.tol=1e-1

info.LambdaSym = [-1, 1, 0]

# info.cont_type = "translation"
# if info.cont_type == "rotation"
#     info.LambdaSym = [-1, 1, 0, 0]
#     info.rotation_center = [0.2, 0.2]
#     info.ncont = 2
#     info.ncont_act = 1 # Two parameters are introduced, but one of them is determined by the other
#     info.params = [0, 0]
# else
#     info.LambdaSym = [-1, 1, 0]
#     info.ncont = 1
#     info.ncont_act = 1
#     info.params = [0]
# end
