
function matcont(nodes::Vector{Snode},M::SparseMatrixCSC{Float64},V::Matrix{ComplexF64},P::Vector{Parametrisation})

println("Init Matcont")

howmany=0
for p in 1:info.max_order
  for i in 1:P[p].m
    howmany+=1
  end
end

Avector=zeros(howmany,info.nza+1)
fdyn=zeros(howmany,info.nza)
fdynpol=zeros(howmany,info.nza)
mappings=zeros(ComplexF64,howmany,info.NN,3)
mappings_r=zeros(howmany,info.NN,3)
u0=zeros(info.NN,3)
mappings_p=zeros(ComplexF64,howmany,info.NN,1)
mappings_pr=zeros(howmany,info.NN,1)
mappings_modal=zeros(howmany,info.neig)
mappings_modal_p=zeros(howmany,info.neig)

index=1
for p in 1:info.max_order
  for i in 1:P[p].m
    Avector[index,:]=P[p].Av[i]
    index+=1
  end
end

index=1
for p in 1:info.max_order
  for i in 1:P[p].m
    for j in 1:info.nmm
      fdyn[index,j]=2*real(P[p].fr[j,i])
    end
    for j in 1:info.nmm
      fdyn[index,j+info.nmm]=-2*imag(P[p].fr[j,i])
    end
    index+=1
  end
end

index=1
# fdynpol[index,j] contains the ̇ρ equation, fdynpol[index,j+info.nmm] contains the ω equation
for p in 1:info.max_order
  for i in 1:P[p].m
    for j in 1:info.nmm
      fdynpol[index,j]=real(P[p].f[j,i])
    end
    for j in 1:info.nmm
      fdynpol[index,j+info.nmm]=imag(P[p].f[j,i])
    end
    index+=1
  end
end

index=1
for p in 1:info.max_order
  for i in 1:P[p].m
    for inode=1:info.NN
      for idof=1:2
        dof=nodes[inode].udof[idof]
        if dof>0
          mappings[index,inode,idof]=P[p].W[dof,i]
          mappings_r[index,inode,idof]=real(P[p].Wr[dof,i])
        end
        u0[inode,idof]=nodes[inode].u[idof]

        dof=nodes[inode].pdof
        if dof>0
          mappings_p[index,inode,1]=P[p].W[dof,i]
          mappings_pr[index,inode,1]=real(P[p].Wr[dof,i])
        end
      end
    end
#    for imode=1:info.neig # all computed modes
#      mappings_modal_vel[index,imode]=V[:,imode]'*M*real.(P[p].Wr[1:info.neq,i])
#      mappings_modal[index,imode]=V[:,imode]'*M*real.(P[p].Wr[info.neq+1:2*info.neq,i])
#    end
    index+=1
  end
end

file = matopen("./output/param.mat","w")
write(file,"M",M)
write(file, "ndof",info.neq)
write(file, "nz",info.nza)
write(file, "max_order",info.max_order)
write(file, "mappings",mappings)
write(file,"mappings_r",mappings_r)
write(file, "u0", u0)
#write(file, "mappings_modal",mappings_modal)
write(file, "mappings_p",mappings_p)
write(file, "mappings_pr",mappings_pr)
#write(file, "mappings_modal_p",mappings_modal_vel)
write(file, "Avector",Avector)
write(file, "fdyn",fdyn)
write(file, "fdynpol",fdynpol)
close(file)

end

