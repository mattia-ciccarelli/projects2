#
# sets input
#
println(Threads.nthreads())
#input_file="demo.jl"
# input_file="cylinder.jl"
# input_file="cylinder_coarse.jl"
# input_file="cylinder_sym.jl"
# input_file="cylinder_adim.jl"
# input_file="cylinder2.jl"
# input_file="fluidic_pinball.jl"
input_file="U_channel_cop.jl"
# input_file="two_cylinders.jl"

using Distributed
using KrylovKit

using SparseArrays
using LinearAlgebra
using Arpack
using ExtendableSparse
const SparseMatrixLNK = ExtendableSparse.SparseMatrixLNK 
using BenchmarkTools
using Pardiso
using MAT
using MATLAB
using Debugger
using UnicodePlots

include("./source/NS2D_defs.jl")
include("./source/NS2D_globalprocedure.jl")
include("./source/NS2D_globalprocedure_iterativeRe.jl")
include("./source/NS2D_readgmsh.jl")
include("./source/NS2D_outgmsh.jl")
include("./source/NS2D_outvtk.jl")
include("./source/NS2D_elemental.jl")
include("./source/NS2D_analysis.jl")
include("./source/NS2D_eig.jl")
include("./source/NS2D_param_struct.jl")
include("./source/NS2D_dpim.jl")
include("./source/NS2D_dpim_FEM.jl")
include("./source/NS2D_dpim_realification.jl")
include("./source/NS2D_dpim_output.jl")
include("./source/quadrature.jl")
include("./source/shape_functions.jl")
include("./source/NS2D_matcont.jl")
include("./source/NS2D_backwardEuler_opt.jl")
include("./source/NS2D_backwardEuler_opt_nou0.jl")
include("./source/NS2D_crankNicolson_opt_nou0.jl")
include("./source/NS2D_inflow.jl")
include("./source/NS2D_test.jl")

global eigaround = 1.7*im
globalprocedure(input_file)
# globalprocedure_iterativeRe()              
# backwardEuler_opt(input_file)
# backwardEuler_opt_nou0(input_file)
# crankNicolson_opt_nou0(input_file)

# include("./source/NS2D_trackeigs.jl")

# track_eigenvalues(4100,4200,20)
