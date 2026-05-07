function globalprocedure(input_file)


fname = "./input/"*input_file
include(fname)
# info, file, Re, Uin, coeffnorm, material, solid, dbc, dbcval, inflowbc = load_input()
info.mesh_file=file

T6,B3,nodes = readgmsh!(info)

##############################
# static solution and matrix assemblage
##############################
global first_run = true
K,M,F0=analysis(nodes,T6,info)
export_u0(nodes, B3, T6, "./output", "")
#test(nodes,T6,info)

##############################
# eigenvalue computation
##############################

@time D,VR,VL,VR_nu=eigAB_Uchannel_parall2(K,M,F0,nodes)

export_eig(nodes, B3, T6, D, real(VR[1:info.neq, :]), "./output","")
#export_eig(nodes, B3, T6, D, real(M*VR[1:info.neq, :]), "./output","")
export_eig2(nodes, B3, T6, D, imag(VR[1:info.uneq, :]), "./output","")
export_eig_nu(nodes, B3, T6, D, real(VR_nu[1:info.uneq, :]), "./output", "_nu")

##############################
# launch DPIM
##############################
P=dpim(K,M,D,VR,VL,VR_nu,nodes,T6,B3) # computes parametrization
realification!(P)  # performs realification

##################################################
#  OUTPUT ON FILE
##################################################

output(P)  # output

##################################################
#  output for matcont
##################################################

matcont(nodes,M,VR,P)

end

