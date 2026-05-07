
function test(nodes::Vector{Snode},T6::Vector{ST6},info::Sinfo)

# assignes inflow
for n in 1:info.NN 
  for d in 1:2       
    dof=nodes[n].udof[d]
    if dof==-2                      
      nodes[n].u[d]=inflow(nodes[n].coor,d)
    end 
  end  
end

uneq=0
# equation numbering: disp 
for e in 1:info.NE           
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      for d in 1:2
        if nodes[n].udof[d]==0   
          uneq+=1             
          nodes[n].udof[d]=uneq   
        end  
      end
    end
  end   
end

pneq=0
# equation numbering: press 
for e in 1:info.NE           
  mat=T6[e].mat
  if mat>0 
    for k in 1:3
      n=T6[e].nodes[k]   
      if nodes[n].pdof==0   
        pneq+=1             
        nodes[n].pdof=uneq+pneq   
      end  
    end
  end   
end
info.uneq=uneq
info.pneq=pneq
neq=uneq+pneq
info.neq=uneq+pneq


# allocation of elements arrays
Xe=zeros(Float64,6,2)
udofe=zeros(Int64,12,1)
pdofe=zeros(Int64,3,1)
Ue=zeros(Float64,6,2)
Pe=zeros(Float64,3,1)
FeU0=zeros(Float64,12)
FeP0=zeros(Float64,3)
FeP1=zeros(Float64,3)
FeU1=zeros(Float64,12)
FeU2=zeros(Float64,12)
FeUT=zeros(Float64,12)
FePT=zeros(Float64,3)
Ke= zeros(Float64,12,12)
De= zeros(Float64,12,3)

e=1

Ue=rand(6,2)
DUe=rand(6,2)

Pe=0*rand(3,1)
DPe=0*rand(3,1)

vUe=zeros(Float64,12,1)
for k=1:6
 vUe[2*k-1:2*k]=DUe[k,:]
end 

#for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      udofe[(k-1)*2+1:k*2]=nodes[n].udof
    end
    for k in 1:3
      n=T6[e].nodes[k]   
      pdofe[k]=nodes[n].pdof
    end
  end
#end  

println("")
println("")
println("")

T6_KeNS!(Ke,De,FeUT,FePT,Xe,Ue+DUe,Pe+DPe,material[mat])
println("FeUT: ",FeUT)
println("")

T6_KeNS!(Ke,De,FeU0,FeP0,Xe,Ue,Pe,material[mat])
println("FeU0: ",FeU0)
println("")

FeU1[:]+=(Ke*vUe+De*DPe)  
FeP1[:]+=De'*vUe
println("FeU1: ",FeU1)
println("")

T6_quad!(FeU2,Xe,DUe)
println("FeU2: ",FeU2)
println("")

println("FeUT2: ",FeU0+FeU1+FeU2)
println("")
println("")


#println("FePT: ",FePT)
#println("")
#println("FePT2: ",FeP0+FeP1)


end   


