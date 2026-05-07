function globalprocedure_iterativeRe()


fname = "./input/"*input_file
include(fname)
info.mesh_file=file

T6,B3,nodes = readgmsh!(info)

##############################
# static solution and matrix assemblage
##############################
global first_run = true
K,M,F0=analysis(nodes,T6,info)
global first_run = false
# for Re = 20.0:1.0:120.0
for Re = 70.0:3.0:120.0
  println("\n\nRe = ", Re, "\n\n")
  global material = [[1.0, coeffnorm*0.1/Re]]
  K,M,F0=analysis(nodes,T6,info)

  ##############################
  # eigenvalue computation
  ##############################

  D,VR,VL,VR_nu=eigABiter(K,M,F0,nodes)
  global eigaround = D[1]
  Re_int = lpad(string(Int(round(Re * 100))), 5, '0')
  if !isdir("./output/_iterativeRe/data_cyl/cyl$(Re_int)/eig")
    mkpath("./output/_iterativeRe/data_cyl/cyl$(Re_int)/eig")
  end
  if Re > 100
    export_eig(nodes, B3, T6, D, real(M*VR[1:info.neq, :]), "./output","")
    l
  end
  export_eigvalues(nodes, B3, T6, D, real(VR[1:info.uneq, :]), "./output/_iterativeRe/data_cyl/cyl$(Re_int)","")
  # export_eig_nu(nodes, B3, T6, D, real(VR_nu[1:info.uneq, :]), "./output/cyl$(Re_int)", "_nu")
  # export_u0(nodes, B3, T6, "./output/cyl$(Re_int)", "")
end

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

