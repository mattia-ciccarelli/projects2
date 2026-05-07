function T6_quad!(FeU::Vector{ComplexF64},Xe::Matrix{Float64},Ue1::Matrix{ComplexF64},Ue2::Matrix{ComplexF64})

fill!(FeU,0.0)

dNda=zeros(Float64,6,2)
N6=zeros(Float64,6,1)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
N=zeros(Float64,12,2)

qr=quadrature_points(Val(:TRI6gp))
for (w,a) in qr

  dNda!(dNda,a,Val(:TRI6n))
  N!(N6,a,Val(:TRI6n))
  Jac[:,:] = Xe'*dNda
  J =  Jac[1,1]*Jac[2,2]-Jac[2,1]*Jac[1,2]
  invJac[:,:] =[Jac[2,2] -Jac[1,2];
             -Jac[2,1] Jac[1,1]]/J
  dNdx[:,:] = dNda*invJac
  for i = 1:6
    ipos1=2*i-1
    ipos2=ipos1+1
    N[ipos1,1] = N6[i]
    N[ipos2,2] = N6[i]
  end
  GU1=Ue1'*dNdx
  GU2=Ue2'*dNdx
  U1=Ue1'*N6
  U2=Ue2'*N6

  coeff=J*w/2
  FeU[:]+=N*[GU1[1,1]*U2[1]+GU1[1,2]*U2[2]; GU1[2,1]*U2[1]+GU1[2,2]*U2[2]]*coeff
  FeU[:]+=N*[GU2[1,1]*U1[1]+GU2[1,2]*U1[2]; GU2[2,1]*U1[1]+GU2[2,2]*U1[2]]*coeff

end  # loop gauss points

return nothing

end



function T6_KeNS!(Ke::Matrix{Float64},De::Matrix{Float64},FeU::Vector{Float64},FeP::Vector{Float64},Fe0::Vector{Float64},
                Xe::Matrix{Float64},Ue::Matrix{Float64},Pe::Matrix{Float64},rho::Float64,nu::Float64,coeffnorm::Float64)

fill!(Ke,0.0)
fill!(De,0.0)
fill!(FeU,0.0)
fill!(FeP,0.0)
fill!(Fe0,0.0)

dNda=zeros(Float64,6,2)
N6=zeros(Float64,6,1)
Jac=zeros(Float64,2,2)
invJac=zeros(Float64,2,2)
dNdx=zeros(Float64,6,2)
B=zeros(Float64,12,3)
N=zeros(Float64,12,2)
G=zeros(Float64,12,2,2)
N3=zeros(Float64,1,3)
KNLe=zeros(Float64,12,12)

D=nu*[2 0 0;
      0 2 0;
      0 0 1];

qr=quadrature_points(Val(:TRI6gp))
for (w,a) in qr

  dNda!(dNda,a,Val(:TRI6n))
  N!(N6,a,Val(:TRI6n))
  N3[:]=a[:]
  Jac[:] = Xe'*dNda
  J =  Jac[1,1]*Jac[2,2]-Jac[2,1]*Jac[1,2]
  invJac[:] =[Jac[2,2] -Jac[1,2];
             -Jac[2,1] Jac[1,1]]/J
  dNdx[:,:] = dNda*invJac
  for i = 1:6
    ipos1=2*i-1
    ipos2=ipos1+1
    B[ipos1,1] = dNdx[i,1]
    B[ipos2,2] = dNdx[i,2]
    B[ipos1,3] = dNdx[i,2]
    B[ipos2,3] = dNdx[i,1]
    G[ipos1,1,1] = dNdx[i,1]
    G[ipos2,2,1] = dNdx[i,1]
    G[ipos1,1,2] = dNdx[i,2]
    G[ipos2,2,2] = dNdx[i,2]
    N[ipos1,1] = N6[i]
    N[ipos2,2] = N6[i]
  end

  GU=Ue'*dNdx
  U=Ue'*N6

  coeff=J*w
  KNLe[:,:]+=N*[GU[1,1]*N[:,1]'+GU[1,2]*N[:,2]'; GU[2,1]*N[:,1]'+GU[2,2]*N[:,2]']*coeff
  KNLe[:,:]+=N*[G[:,1,1]'*U[1]+G[:,1,2]'*U[2]; G[:,2,1]'*U[1]+G[:,2,2]'*U[2]]*coeff
  FeU[:]+=N*[GU[1,1]*U[1]+GU[1,2]*U[2]; GU[2,1]*U[1]+GU[2,2]*U[2]]*coeff

  coeff=J*w
  Ke[:,:]+=(B*D*B')*coeff/coeffnorm
  coeff=J*w/rho
  De[:,:]-=(B[:,1]+B[:,2])*N3*coeff

end  # loop gauss points

vUe=zeros(Float64,12,1)
for k=1:6
 vUe[2*k-1:2*k]=Ue[k,:]
end

FeU[:]+=(Ke*vUe+De*Pe)
FeP[:]+=De'*vUe
Fe0[:]+=Ke*vUe/nu
Ke[:,:]+=KNLe[:,:]

return nothing

end




# mass matrix: T6
function T6_Me!(Me::Matrix{Float64},X::Matrix{Float64},mate::Vector{Float64})

fill!(Me,0.0)
dNda=zeros(Float64,(6,2))
NL=zeros(Float64,6)

rho=mate[1]
qr=quadrature_points(Val(:TRI6gp))

for (w,a) in qr
  dNda!(dNda,a,Val(:TRI6n))
  F=X'*dNda
  J=F[1,1]*F[2,2]-F[1,2]*F[2,1]
  N!(NL,a,Val(:TRI6n))
  N=[NL[1] 0 NL[2] 0 NL[3] 0 NL[4] 0 NL[5] 0 NL[6] 0;
     0 NL[1] 0 NL[2] 0 NL[3] 0 NL[4] 0 NL[5] 0 NL[6]]
  Me[:,:]+=(N'*N)*rho*J*w
end

end


