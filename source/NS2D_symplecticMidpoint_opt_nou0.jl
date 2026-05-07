
function symplecticMidpoint_opt_nou0(input_file)

  output_dir = "./output/_timeint/"

  fname = "../DPIM2D_NS/input/" * input_file
  include(fname)
  info.mesh_file = file

  T6, B3, nodes = readgmsh!(info)

  ##############################
  # static solution and matrix assemblage
  ##############################
  global first_run = true
  K, M, F0 = analysis(nodes, T6, info)

  # Correct mass with lifting terms
  Fmass = zeros(Float64, info.neq)
  # M = mass_assembly(nodes, T6, Fmass)

  # Perturb u0 to speed up transient
  mfile = matopen("./output/setup_u0_timeint.mat")
  delta_u0 = read(mfile, "delta_u0")
  delta_p0 = read(mfile, "delta_p0")
  close(mfile)

  neq = info.neq
  uneq = info.uneq
  u_n = zeros(Float64, neq)

  # assignes inflow
  for n in 1:info.NN
    for d in 1:2
      dof = nodes[n].udof[d]
      if dof == -2
        nodes[n].u[d] = inflow(nodes[n].coor, d)
      end
    end
  end

  # Perturb solution
  for n in 1:info.NN
    for d in 1:2
      udof = nodes[n].udof[d]
      if udof > 0
        u_n[udof] += delta_u0[n, d]
      end
    end
    pdof = nodes[n].pdof
    if pdof > 0
      u_n[pdof] += delta_p0[n]
    end
  end


  nnu = 1
  nu = material[1][2]
  dnu = -coeffnorm * 0.0001
  vnu = zeros(Float64, nnu)

  omega = 16.85
  period = 2 * pi / omega
  nperiods = 40
  nstepd = 100
  te = nperiods * period
  dt = period / nstepd
  char_dt = 0.1 / Re
  println("dt: ", dt)
  println("Characteristic time: ", char_dt)
  dt = min(dt, 0.7*char_dt)
  nstepd = floor(Int, period / dt)
  nstep = floor(Int, te / dt) + 1

  Mat = SparseMatrixLNK(neq, neq)
  Mat .= 4 * M / dt + 2 * K
  Mat = SparseMatrixCSC(Mat)
  Mat_NR = SparseMatrixLNK(neq, neq)
  Mat_NR = SparseMatrixCSC(Mat_NR)

  Rhs = zeros(Float64, neq)
  Rhs_NR = zeros(Float64, neq)
  u_nplus1 = zeros(Float64, neq)
  u_nplus1 .= u_n
  Du = zeros(Float64, neq)

  for inu = 1:nnu
    re_value = coeffnorm * 10.0 / nu
    timestep_str = lpad(string(0), 5, "0")
    export_timestep(nodes, B3, T6, "$(output_dir)", "/Re$(re_value)_$(timestep_str)", u_n)

    t = 0
    #fill!(u,0.0)

    vnu[inu] = nu

    for i = 1:nstep

      t += dt
      print("\ri_nu: ", inu, ", Period n°: ", ceil(Int, i / nstepd), " of ", nperiods, ", Step ", (i - 1) % nstepd + 1, " of ", nstepd, " | Residuum: N/A                       ")

      residuum = 100
      Rhs .= 4 * M / dt * u_n - 2 * K * u_n
      assembly_RHS_quad_NS!(Rhs, u_n, u_n, nodes, T6)


      while residuum > 2e-8

        Mat_NR .= Mat
        assembly_Clin!(Mat_NR, nodes, T6, u_nplus1)
        assembly_Clin!(Mat_NR, nodes, T6, u_n)

        Rhs_NR .= Rhs - Mat * u_nplus1
        assembly_RHS_quad_NS!(Rhs_NR, u_nplus1, u_nplus1, nodes, T6)
        assembly_RHS_quad_NS!(Rhs_NR, u_n     , u_nplus1, nodes, T6)
        assembly_RHS_quad_NS!(Rhs_NR, u_nplus1, u_n     , nodes, T6)

        # Solve for Du
        Du = Mat_NR \ Rhs_NR
        # Update delu_nplus1
        u_nplus1 .+= Du

        # If norm(Du) < tol, break
        residuum = norm(Du)
        print("\ri_nu: ", inu, ", Period n°: ", ceil(Int, i / nstepd), " of ", nperiods, ", Step ", (i - 1) % nstepd + 1, " of ", nstepd, " | Residuum: ", residuum, "      ")
        # println(residuum)
      end

      # Update u_nplus1
      u_n .= u_nplus1
      # Save u_nplus1, 30 times per period
      if (i % floor(Int, nstepd / 30) == 0)
        timestep_str = lpad(string(i), 5, "0")
        export_timestep(nodes, B3, T6, "$(output_dir)", "/Re$(re_value)_$(timestep_str)", u_n)
      end

    end  # end iterations

    nu += dnu



  end # end iomega

end




function assembly_RHS_quad_NS!(Rhs::Vector{Float64}, U1::Vector{Float64}, U2::Vector{Float64},
  nodes::Vector{Snode}, T6::Vector{ST6})

  Xe = zeros(Float64, 6, 2)
  dofe = zeros(Int64, 12)
  Fe = zeros(ComplexF64, 12)
  U1e = zeros(ComplexF64, 6, 2)
  U2e = zeros(ComplexF64, 6, 2)

  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:6
        n = T6[e].nodes[k]
        Xe[k, :] = nodes[n].coor
        dofe[2*k-1:2*k] = nodes[n].udof
        for idof = 1:2
          dof = nodes[n].udof[idof]
          if dof > 0
            U1e[k, idof] = U1[dof]
            U2e[k, idof] = U2[dof]
          else
            U1e[k, idof] = 0.0
            U2e[k, idof] = 0.0
          end
        end
      end

      T6_RHS_quad_NS!(Fe, Xe, U1e, U2e)

      for i in 1:12
        dofU = dofe[i]
        if dofU > 0
          Rhs[dofU] -= Fe[i]
        end
      end
    end

  end  # loop over NE
end

function T6_RHS_quad_NS!(Fe::Vector{ComplexF64}, Xe::Matrix{Float64},
  U1e::Matrix{ComplexF64}, U2e::Matrix{ComplexF64})

  fill!(Fe, 0.0)

  dNda = zeros(Float64, 6, 2)
  N6 = zeros(Float64, 6, 1)
  Jac = zeros(Float64, 2, 2)
  invJac = zeros(Float64, 2, 2)
  dNdx = zeros(Float64, 6, 2)
  N = zeros(Float64, 12, 2)
  B = zeros(Float64, 12, 3)

  qr = quadrature_points(Val(:TRI6gp))
  for (w, a) in qr

    dNda!(dNda, a, Val(:TRI6n))
    N!(N6, a, Val(:TRI6n))
    Jac[:, :] = Xe' * dNda
    J = Jac[1, 1] * Jac[2, 2] - Jac[2, 1] * Jac[1, 2]
    invJac[:, :] = [
      Jac[2, 2] -Jac[1, 2];
      -Jac[2, 1] Jac[1, 1]] / J
    dNdx[:, :] = dNda * invJac
    for i = 1:6
      ipos1 = 2 * i - 1
      ipos2 = ipos1 + 1
      B[ipos1, 1] = dNdx[i, 1]
      B[ipos2, 2] = dNdx[i, 2]
      B[ipos1, 3] = dNdx[i, 2]
      B[ipos2, 3] = dNdx[i, 1]
      N[ipos1, 1] = N6[i]
      N[ipos2, 2] = N6[i]
    end
    GU1 = transpose(U1e) * dNdx
    GU2 = transpose(U2e) * dNdx
    U1 = transpose(U1e) * N6
    U2 = transpose(U2e) * N6

    coeff = J * w

    Fe[:] += 0.5 * N * ([GU1[1, 1] * U2[1] + GU1[1, 2] * U2[2]; GU1[2, 1] * U2[1] + GU1[2, 2] * U2[2]] +
                        [GU2[1, 1] * U1[1] + GU2[1, 2] * U1[2]; GU2[2, 1] * U1[1] + GU2[2, 2] * U1[2]]) * coeff

  end  # loop gauss points

  return nothing

end

function assembly_Clin!(K::SparseMatrixCSC{Float64}, nodes::Vector{Snode}, T6::Vector{ST6}, u::Vector{Float64})

  # allocation of elements arrays
  Xe = zeros(Float64, 6, 2)
  udofe = zeros(Int64, 12, 1)
  Ue = zeros(Float64, 6, 2)
  Ce = zeros(Float64, 12, 12)

  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:6
        n = T6[e].nodes[k]
        Xe[k, :] = nodes[n].coor
        udofe[(k-1)*2+1:k*2] = nodes[n].udof

        for idim in 1:2
          dofi = nodes[n].udof[idim]
          if dofi > 0
            Ue[k, idim] = u[dofi]
          else
            Ue[k, idim] = 0.0
          end
        end
      end
      T6_CeLin!(Ce, Xe, Ue)
      for i = 1:12
        dofi = udofe[i]
        if dofi > 0
          for j = 1:12
            dofj = udofe[j]
            if dofj > 0
              K[dofi, dofj] += Ce[i, j]
            end
          end
        end
      end
    end
  end
end


function T6_CeLin!(Ce::Matrix{Float64}, Xe::Matrix{Float64}, Ue::Matrix{Float64})
  fill!(Ce, 0.0)

  dNda = zeros(Float64, 6, 2)
  N6 = zeros(Float64, 6, 1)
  Jac = zeros(Float64, 2, 2)
  invJac = zeros(Float64, 2, 2)
  dNdx = zeros(Float64, 6, 2)
  B = zeros(Float64, 12, 3)
  N = zeros(Float64, 12, 2)
  G = zeros(Float64, 12, 2, 2)

  qr = quadrature_points(Val(:TRI6gp))
  for (w, a) in qr

    dNda!(dNda, a, Val(:TRI6n))
    N!(N6, a, Val(:TRI6n))
    Jac[:] = Xe' * dNda
    J = Jac[1, 1] * Jac[2, 2] - Jac[2, 1] * Jac[1, 2]
    invJac[:] = [
      Jac[2, 2] -Jac[1, 2];
      -Jac[2, 1] Jac[1, 1]] / J
    dNdx[:, :] = dNda * invJac
    for i = 1:6
      ipos1 = 2 * i - 1
      ipos2 = ipos1 + 1
      B[ipos1, 1] = dNdx[i, 1]
      B[ipos2, 2] = dNdx[i, 2]
      B[ipos1, 3] = dNdx[i, 2]
      B[ipos2, 3] = dNdx[i, 1]
      G[ipos1, 1, 1] = dNdx[i, 1]
      G[ipos2, 2, 1] = dNdx[i, 1]
      G[ipos1, 1, 2] = dNdx[i, 2]
      G[ipos2, 2, 2] = dNdx[i, 2]
      N[ipos1, 1] = N6[i]
      N[ipos2, 2] = N6[i]
    end

    GU = Ue' * dNdx
    U = Ue' * N6

    coeff = J * w
    Ce[:, :] += N * [GU[1, 1] * N[:, 1]' + GU[1, 2] * N[:, 2]'; GU[2, 1] * N[:, 1]' + GU[2, 2] * N[:, 2]'] * coeff
    Ce[:, :] += N * [G[:, 1, 1]' * U[1] + G[:, 1, 2]' * U[2]; G[:, 2, 1]' * U[1] + G[:, 2, 2]' * U[2]] * coeff

  end  # loop gauss points

  return nothing
end



function mass_assembly(nodes::Vector{Snode}, T6::Vector{ST6}, Fmass::Vector{Float64})

  neq = info.neq
  M = SparseMatrixLNK(neq, neq)

  # allocation of elements arrays
  Xe = zeros(Float64, 6, 2)
  dofe = zeros(Int64, 12, 1)
  Me = zeros(Float64, 12, 12)

  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:6
        n = T6[e].nodes[k]
        Xe[k, :] = nodes[n].coor
        dofe[(k-1)*2+1:k*2] = nodes[n].udof
      end
      T6_Me!(Me, Xe, material[mat])
      for i = 1:12
        dofi = dofe[i]
        if dofi > 0
          for j = 1:12
            dofj = dofe[j]
            if dofj > 0
              M[dofi, dofj] += Me[i, j]
            end
          end
        end
      end
    end
  end

  M = SparseMatrixCSC(M)
  return M
end



function export_timestep(Lnode::Vector{Snode}, Lel1D::Vector{SB3}, Sel2D::Vector{ST6}, out_dir::String, tag::String, u)
  odir = out_dir
  #
  try
    mkdir(odir)
  catch
  end
  #
  oeigf = odir * tag * ".vtk"
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
  write(io, "VECTORS u float\n")
  for i = 1:info.NN
    for j = 1:2
      if Lnode[i].udof[j] > 0
        write(io, string(u[Lnode[i].udof[j]]) * " ")
      else
        write(io, string(0.0) * " ")
      end
    end
    write(io, string(0.0))
    write(io, '\n')
  end
  write(io, '\n')
  write(io, "SCALARS p float\n")
  write(io, "LOOKUP_TABLE default\n")
  for i = 1:info.NN
    if Lnode[i].pdof > 0
      write(io, string(u[Lnode[i].pdof]) * " ")
    else
      write(io, string(0.0) * " ")
    end
    write(io, '\n')
  end
  #
  close(io)
  #
  return nothing
end
