function export_eigvalues(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, D::Vector{ComplexF64}, VR::Matrix{Float64},
  out_dir::String, tag::String)
  odir = out_dir * "/eig"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigv = odir * "/eigenfrequencies" * tag * ".txt"
  #
  io = open(oeigv, "w")
  #
  for i = 1:size(D)[1]
    write(io, string(D[i]))
    write(io, "\n")
  end
  close(io)
  #
  return nothing
end

function export_eig(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, D::Vector{ComplexF64}, VR::Matrix{Float64},
  out_dir::String, tag::String)
  odir = out_dir * "/eig"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigv = odir * "/eigenfrequencies" * tag * ".txt"
  oeigf = odir * "/eigenfunctions" * tag * ".vtk"
  #
  io = open(oeigv, "w")
  #
  for i = 1:size(D)[1]
    write(io, string(D[i]))
    write(io, "\n")
  end
  close(io)
  #
  io = open(oeigf, "w")
  write(io, "# vtk DataFile Version 3.0\n")
  write(io, "eigenfunctions\n")
  write(io, "ASCII\n")
  write(io, "\n")
  write(io, "DATASET UNSTRUCTURED_GRID\n")
  write(io, "POINTS " * string(info.NN) * " float\n")
  #
  for n = 1:info.NN
    for i = 1:2
      write(io, string(Lnode[n].coor[i]) * " ")
    end
    write(io, string(0.0))
    write(io, "\n")
  end
  write(io, "\n")
  # compute total cell number and entries
  tcn = info.NE
  tce = tcn * 6
  #
  write(io, "CELLS " * string(tcn) * " " * string(tce + tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "6 ")
    write(io, string(Sel2D[e].nodes[1] - 1) * " ")
    write(io, string(Sel2D[e].nodes[2] - 1) * " ")
    write(io, string(Sel2D[e].nodes[3] - 1) * " ")
    write(io, string(Sel2D[e].nodes[4] - 1) * " ")
    write(io, string(Sel2D[e].nodes[5] - 1) * " ")
    write(io, string(Sel2D[e].nodes[6] - 1))
    write(io, "\n")
  end
  #
  write(io, "\n")
  write(io, "CELL_TYPES " * string(tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "22\n")
  end
  #
  write(io, "\n")
  write(io, "POINT_DATA " * string(info.NN) * "\n")
  #
  for eigm = 1:size(D)[1]
    write(io, "VECTORS eigenmode" * lpad(string(eigm), 3, '0') * " float\n")
    for i = 1:info.NN
      for j = 1:2
        if Lnode[i].udof[j] > 0
          write(io, string(VR[Lnode[i].udof[j], eigm]) * " ")
        else
          write(io, string(0.0) * " ")
        end
      end
      write(io, string(0.0))
      write(io, '\n')
    end
    write(io, '\n')
  end
  #
  close(io)
  #
  return nothing
end

function export_eig_nu(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, D::Vector{ComplexF64}, VR_nu::Matrix{Float64},
  out_dir::String, tag::String)
  odir = out_dir * "/eig"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigv = odir * "/eigenfrequencies" * tag * ".txt"
  oeigf = odir * "/eigenfunctions" * tag * ".vtk"
  #
  io = open(oeigv, "w")
  #
  for i = 1:size(D)[1]
    write(io, string(D[i]))
    write(io, "\n")
  end
  close(io)
  #
  io = open(oeigf, "w")
  write(io, "# vtk DataFile Version 3.0\n")
  write(io, "eigenfunctions\n")
  write(io, "ASCII\n")
  write(io, "\n")
  write(io, "DATASET UNSTRUCTURED_GRID\n")
  write(io, "POINTS " * string(info.NN) * " float\n")
  #
  for n = 1:info.NN
    for i = 1:2
      write(io, string(Lnode[n].coor[i]) * " ")
    end
    write(io, string(0.0))
    write(io, "\n")
  end
  write(io, "\n")
  # compute total cell number and entries
  tcn = info.NE
  tce = tcn * 6
  #
  write(io, "CELLS " * string(tcn) * " " * string(tce + tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "6 ")
    write(io, string(Sel2D[e].nodes[1] - 1) * " ")
    write(io, string(Sel2D[e].nodes[2] - 1) * " ")
    write(io, string(Sel2D[e].nodes[3] - 1) * " ")
    write(io, string(Sel2D[e].nodes[4] - 1) * " ")
    write(io, string(Sel2D[e].nodes[5] - 1) * " ")
    write(io, string(Sel2D[e].nodes[6] - 1))
    write(io, "\n")
  end
  #
  write(io, "\n")
  write(io, "CELL_TYPES " * string(tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "22\n")
  end
  #
  write(io, "\n")
  write(io, "POINT_DATA " * string(info.NN) * "\n")
  #
  write(io, "VECTORS eigenmode_nu float\n")
  for i = 1:info.NN
    for j = 1:2
      if Lnode[i].udof[j] > 0
        write(io, string(VR_nu[Lnode[i].udof[j]]) * " ")
      else
        write(io, string(0.0) * " ")
      end
    end
    write(io, string(0.0))
    write(io, '\n')
  end
  write(io, '\n')
  #
  close(io)
  #
  return nothing
end




function print_cell_nodes(io, iΩ, set, nn, skip, ::Type{Val{:P18}})
  #
  for e = 1:iΩ.ne[set]
    write(io, "18 ")
    write(io, string(iΩ.e2n[skip+1+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+2+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+3+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+4+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+5+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+6+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+7+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+8+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+9+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+10+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+11+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+12+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+13+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+14+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+15+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+16+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+17+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+18+(e-1)*nn] - 1))
    write(io, "\n")
  end
  #
  return nothing
end





function export_u0(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, out_dir::String, tag::String)
  odir = out_dir * "/eig"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigf = odir * "/uu0" * tag * ".vtk"
  #
  io = open(oeigf, "w")
  write(io, "# vtk DataFile Version 3.0\n")
  write(io, "u0\n")
  write(io, "ASCII\n")
  write(io, "\n")
  write(io, "DATASET UNSTRUCTURED_GRID\n")
  write(io, "POINTS " * string(info.NN) * " float\n")
  #
  for n = 1:info.NN
    for i = 1:2
      write(io, string(Lnode[n].coor[i]) * " ")
    end
    write(io, string(0.0))
    write(io, "\n")
  end
  write(io, "\n")
  # compute total cell number and entries
  tcn = info.NE
  tce = tcn * 6
  #
  write(io, "CELLS " * string(tcn) * " " * string(tce + tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "6 ")
    write(io, string(Sel2D[e].nodes[1] - 1) * " ")
    write(io, string(Sel2D[e].nodes[2] - 1) * " ")
    write(io, string(Sel2D[e].nodes[3] - 1) * " ")
    write(io, string(Sel2D[e].nodes[4] - 1) * " ")
    write(io, string(Sel2D[e].nodes[5] - 1) * " ")
    write(io, string(Sel2D[e].nodes[6] - 1))
    write(io, "\n")
  end
  #
  write(io, "\n")
  write(io, "CELL_TYPES " * string(tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "22\n")
  end
  #
  write(io, "\n")
  write(io, "POINT_DATA " * string(info.NN) * "\n")
  #
  write(io, "VECTORS u0 float\n")
  for i = 1:info.NN
    for j = 1:2
      if Lnode[i].udof[j] > 0
        write(io, string(Lnode[i].u[j]) * " ")
      else
        write(io, string(Lnode[i].u[j]) * " ")
      end
    end
    write(io, string(0.0))
    write(io, '\n')
  end
  write(io, '\n')
  #
  close(io)
  #
  return nothing
end



function export_du(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, out_dir::String, tag::String, du)
  odir = out_dir * "/eig"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigf = odir * "/du" * tag * ".vtk"
  #
  io = open(oeigf, "w")
  write(io, "# vtk DataFile Version 3.0\n")
  write(io, "u0\n")
  write(io, "ASCII\n")
  write(io, "\n")
  write(io, "DATASET UNSTRUCTURED_GRID\n")
  write(io, "POINTS " * string(info.NN) * " float\n")
  #
  for n = 1:info.NN
    for i = 1:2
      write(io, string(Lnode[n].coor[i]) * " ")
    end
    write(io, string(0.0))
    write(io, "\n")
  end
  write(io, "\n")
  # compute total cell number and entries
  tcn = info.NE
  tce = tcn * 6
  #
  write(io, "CELLS " * string(tcn) * " " * string(tce + tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "6 ")
    write(io, string(Sel2D[e].nodes[1] - 1) * " ")
    write(io, string(Sel2D[e].nodes[2] - 1) * " ")
    write(io, string(Sel2D[e].nodes[3] - 1) * " ")
    write(io, string(Sel2D[e].nodes[4] - 1) * " ")
    write(io, string(Sel2D[e].nodes[5] - 1) * " ")
    write(io, string(Sel2D[e].nodes[6] - 1))
    write(io, "\n")
  end
  #
  write(io, "\n")
  write(io, "CELL_TYPES " * string(tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "22\n")
  end
  #
  write(io, "\n")
  write(io, "POINT_DATA " * string(info.NN) * "\n")
  #
  write(io, "VECTORS du float\n")
  for i = 1:info.NN
    for j = 1:2
      if Lnode[i].udof[j] > 0
        write(io, string(du[Lnode[i].udof[j]]) * " ")
      else
        write(io, string(0.0) * " ")
      end
    end
    write(io, string(0.0))
    write(io, '\n')
  end
  write(io, '\n')
  #
  close(io)
  #
  return nothing
end




function print_cell_nodes(io, iΩ, set, nn, skip, ::Type{Val{:P18}})
  #
  for e = 1:iΩ.ne[set]
    write(io, "18 ")
    write(io, string(iΩ.e2n[skip+1+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+2+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+3+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+4+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+5+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+6+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+7+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+8+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+9+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+10+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+11+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+12+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+13+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+14+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+15+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+16+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+17+(e-1)*nn] - 1) * " ")
    write(io, string(iΩ.e2n[skip+18+(e-1)*nn] - 1))
    write(io, "\n")
  end
  #
  return nothing
end



function export_eig2(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, D::Vector{ComplexF64}, VR::Matrix{Float64},
  out_dir::String, tag::String)
  odir = out_dir * "/eig2"
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigv = odir * "/eigenfrequencies" * tag * ".txt"
  oeigf = odir * "/eigenfunctions" * tag * ".vtk"
  #
  io = open(oeigv, "w")
  #
  for i = 1:size(D)[1]
    write(io, string(D[i]))
    write(io, "\n")
  end
  close(io)
  #
  io = open(oeigf, "w")
  write(io, "# vtk DataFile Version 3.0\n")
  write(io, "eigenfunctions\n")
  write(io, "ASCII\n")
  write(io, "\n")
  write(io, "DATASET UNSTRUCTURED_GRID\n")
  write(io, "POINTS " * string(info.NN) * " float\n")
  #
  for n = 1:info.NN
    for i = 1:2
      write(io, string(Lnode[n].coor[i]) * " ")
    end
    write(io, string(0.0))
    write(io, "\n")
  end
  write(io, "\n")
  # compute total cell number and entries
  tcn = info.NE
  tce = tcn * 6
  #
  write(io, "CELLS " * string(tcn) * " " * string(tce + tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "6 ")
    write(io, string(Sel2D[e].nodes[1] - 1) * " ")
    write(io, string(Sel2D[e].nodes[2] - 1) * " ")
    write(io, string(Sel2D[e].nodes[3] - 1) * " ")
    write(io, string(Sel2D[e].nodes[4] - 1) * " ")
    write(io, string(Sel2D[e].nodes[5] - 1) * " ")
    write(io, string(Sel2D[e].nodes[6] - 1))
    write(io, "\n")
  end
  #
  write(io, "\n")
  write(io, "CELL_TYPES " * string(tcn) * "\n")
  #
  for e = 1:info.NE
    write(io, "22\n")
  end
  #
  write(io, "\n")
  write(io, "POINT_DATA " * string(info.NN) * "\n")
  #
  for eigm = 1:size(D)[1]
    write(io, "VECTORS eigenmode" * lpad(string(eigm), 3, '0') * " float\n")
    for i = 1:info.NN
      for j = 1:2
        if Lnode[i].udof[j] > 0
          write(io, string(VR[Lnode[i].udof[j], eigm]) * " ")
        else
          write(io, string(0.0) * " ")
        end
      end
      write(io, string(0.0))
      write(io, '\n')
    end
    write(io, '\n')
  end
  #
  close(io)
  #
  return nothing
end