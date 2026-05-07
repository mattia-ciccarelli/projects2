
function readgmsh!(info::Sinfo)

file=info.mesh_file
fname = "./input/"*file

 # extract mesh infos
fhand=open(fname,"r")

line = readline(fhand)
line = readline(fhand)
line = readline(fhand)
line = readline(fhand)

# reads node number and coordinates and saves them in
# temp variables nodenum and coor
line = readline(fhand)
NN = Meta.parse(line)

coor=zeros(Float64,NN,2)
nodenum=zeros(Int64,NN)
for i in 1:NN
 line=readline(fhand)
 line=split(line)
 nodenum[i]=Meta.parse(line[1])
 coor[i,1]=Meta.parse(line[2])
 coor[i,2]=Meta.parse(line[3])
# coor[i,2]=Meta.parse(line[2])   # pi/2 ccw rotation
# coor[i,1]=-Meta.parse(line[3])
end

maxnodenum=maximum(nodenum);           # max node number in GMSH numbering

nodes= Vector{Snode}(undef,maxnodenum)
for i = 1:NN
  n=nodenum[i];
  nodes[n] = Snode(coor[i,:],zeros(Int64,2),zeros(Float64,2),0,0.0)
end

line = readline(fhand)
line = readline(fhand)

lendbc=0
if (@isdefined dbc)
 lendbc=length(dbc)
end
leninflowbc=0
if (@isdefined inflowbc)
 leninflowbc=length(inflowbc)
end
lentbc=0
if(@isdefined tbc)
 lentbc=length(tbc)
end
lensolid=length(solid)

line = readline(fhand)     # $Elements
# reads elements and allocates list of elements
NET = Meta.parse(line)
listel=zeros(Int64,NET,9)
NE=0
NL=0

for i in 1:NET                 # loop over lines of section $Nodes
#  global NE,NL                    # integers are not global; vectors are. in for loops all var all local
#  local line
  line=readline(fhand)
  line=split(line)
  itemp=Meta.parse(line[2]);          # type of element
  listel[i,2]=itemp
  listel[i,3]=Meta.parse(line[4]);          # physical set
  if itemp == 8
    for j in 4:6
      listel[i,j]=Meta.parse(line[j+2])
    end
    for j in 1:lentbc
#      if tbc[j][1]==listel[i,3]
      if tbc[j]==listel[i,3]
        NL+=1
        listel[i,1]=NL                  # saves the counter in first place
        break
      end
    end
  elseif itemp==9                      # if element is a T6
    for j in 4:9
      listel[i,j]=Meta.parse(line[j+2])
    end
    NE+=1
    listel[i,1]=NE                   # saves the counter in first place
  else
    println("wrong element type: ",itemp)
    return
  end
end

#analysis=Sanalysis(NN,NE,NL,0)
info.NN=NN
info.NE=NE
info.NL=NL

T6= Vector{ST6}(undef,NE)
for i = 1:NE
 T6[i] = ST6(0,zeros(Int64,6))
end

B3= Vector{SB3}(undef,NL)
for i = 1:NL
 B3[i] = SB3(zeros(Int64,3),zeros(Float64,2),0)
end

for i in 1:NET

  iE=listel[i,1]
  if listel[i,2]==8 && iE>0             # B3 loaded element

   physet=listel[i,3]                   # gets physet
   for j in 1:lentbc
#    if tbc[j][1]==physet
    if tbc[j]==physet
      B3[iE].nodes=listel[i,4:6]
#      B3[iE].press=tbc[j][2]
       B3[iE].press=tbcval[j]
    end
   end

  elseif listel[i,2]==8                 # if a dbc element

    physet=listel[i,3]                   # gets physet

    for j in 1:leninflowbc
      if inflowbc[j][1]==physet
        for k in 4:6
          n=listel[i,k]
          compnt=dbc[j][2]
          nodes[n].udof[compnt]=-2
        end
      end
    end

    for j in 1:lendbc
      if dbc[j][1]==physet
        for k in 4:6
          n=listel[i,k]
          compnt=dbc[j][2]
          nodes[n].udof[compnt]=-1
          nodes[n].u[compnt]=dbcval[j]
        end
      end
    end

  elseif listel[i,2]==9                 # if element is a T6

   T6[iE].nodes=listel[i,4:9]
   physet=listel[i,3]

   for j in 1:lensolid
    if solid[j][1]==physet
      T6[iE].mat=solid[j][2]
    end
   end

  end

end

if (@isdefined dbcn)
  ndbcn=length(dbcn)
  for i in 1:ndbcn
   nodes[dbcn[i][1]].dof[dbcn[i][2]]=-1
   nodes[dbcn[i][1]].u[dbcn[i][2]]=dbcnval[i]
  end
end


return T6,B3,nodes
end # end function

