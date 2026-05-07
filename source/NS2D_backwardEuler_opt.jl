
function backwardEuler_opt(input_file)

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
  M = mass_assembly(nodes, T6, Fmass)

  # Perturb u0 to speed up transient
  mfile = matopen("./output/setup_u0_timeint.mat")
  delta_u0 = read(mfile, "delta_u0")
  delta_p0 = read(mfile, "delta_p0")
  close(mfile)

  for n in 1:info.NN
    for d in 1:2
      dof = nodes[n].udof[d]
      nodes[n].u[d] += delta_u0[n, d]
    end
    dof = nodes[n].pdof
    nodes[n].p += delta_p0[n]
  end

  neq = info.neq
  uneq = info.uneq
  u_n = zeros(Float64, neq, 1)

  # assignes inflow
  for n in 1:info.NN
    for d in 1:2
      dof = nodes[n].udof[d]
      if dof == -2
        nodes[n].u[d] = inflow(nodes[n].coor, d)
      end
    end
  end

  # fill dofs
  for n in 1:info.NN
    for d in 1:2
      dofU = nodes[n].udof[d]
      if dofU > 0
        u_n[dofU] = nodes[n].u[d]
      end
    end
    dofP = nodes[n].pdof
    if dofP > 0
      u_n[dofP] = nodes[n].p
    end
  end


  dofnode = nodes[2055].udof[2]

  # Internal forces
  F = zeros(Float64, neq, 1)

  nnu = 1
  nu = material[1][2]
  dnu = -coeffnorm * 0.0001
  vmax = zeros(Float64, nnu)
  vmin = zeros(Float64, nnu)
  vnu = zeros(Float64, nnu)

  omega = 16.85
  period = 2 * pi / omega
  nperiods = 80
  nstepd = 200
  te = nperiods * period

  Rhs = zeros(Float64, neq, 1)
  du = zeros(Float64, neq, 1)
  u_nplus1 = zeros(Float64, neq, 1)
  u_nplus1 .= u_n

  for inu = 1:nnu
    re_value = coeffnorm * 1.0 / nu
    timestep_str = lpad(string(0), 5, "0")
    export_timestep(nodes, B3, T6, "$(output_dir)", "/Re$(re_value)_$(timestep_str)", u_n)

    t = 0
    #fill!(u,0.0)

    vnu[inu] = nu

    dt = period / nstepd
    nstep = floor(Int, te / dt) + 1
    hist = zeros(Float64, nstep, 1)
    for i = 1:nstep

      t += dt
      print("\ri_nu: ", inu, ", Period n°: ", ceil(Int, i / nstepd), " of ", nperiods, ", Step ", (i - 1) % nstepd + 1, " of ", nstepd, " | Residuum: N/A                       ")

      residuum = 100
      iter = 0


      # Reassignes inflow
      #inflow_scale = t < 0 ? 0.5 * (1 - cos(pi * t)) : 1
      inflow_scale = 1
      for n in 1:info.NN
        for d in 1:2
          dof = nodes[n].udof[d]
          if dof == -2
            nodes[n].u[d] = inflow(nodes[n].coor, d) * inflow_scale
          end
        end
      end

      Rhs[:] = M * u_n

      while residuum > 2e-8

        iter += 1
        K = SparseMatrixLNK(neq, neq)
        compute_residual_opt!(F, K, nodes, T6, B3, nu, u_nplus1)

        # Solution phase
        K = SparseMatrixCSC(K)  # convert to other sparse format

        du[:] = (M + dt * K) \ (Rhs - M * u_nplus1 + dt * F)

        u_nplus1 .+= du

        residuum = norm(du)
        print("\ri_nu: ", inu, ", Period n°: ", ceil(Int, i / nstepd), " of ", nperiods, ", Step ", (i - 1) % nstepd + 1, " of ", nstepd, " | Residuum: ", residuum, "      ")
      end

      u_n .= u_nplus1

      hist[i] = u_n[dofnode]
      #export_u0(nodes, B3, T6, "./output", "_$i")
      # export_du(nodes, B3, T6, "./output", "$iter", u_n)

      # Save the entire field `u_n` as a VTK file for this timestep
      timestep_str = lpad(string(i), 5, "0")
      export_timestep(nodes, B3, T6, "$(output_dir)", "/Re$(re_value)_$(timestep_str)", u_n)


    end  # end iterations

    nu += dnu
    println(nu)



  end # end iomega

end


function compute_residual_opt!(F::Matrix{Float64}, K::SparseMatrixLNK{Float64},
  nodes::Vector{Snode}, T6::Vector{ST6}, B3::Vector{SB3}, nu::Float64, U::Matrix{Float64})

  fill!(F, 0.0)

  # allocation of elements arrays
  Xe = zeros(Float64, 6, 2)
  udofe = zeros(Int64, 12, 1)
  pdofe = zeros(Int64, 3, 1)
  Ue = zeros(Float64, 6, 2)
  Pe = zeros(Float64, 3, 1)
  FeU = zeros(Float64, 12)
  Fe0 = zeros(Float64, 12)
  FeP = zeros(Float64, 3)
  Ke = zeros(Float64, 12, 12)
  De = zeros(Float64, 12, 3)

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
            Ue[k, idim] = U[dofi]
          else
            Ue[k, idim] = nodes[n].u[idim]
          end
        end
      end
      for k in 1:3
        n = T6[e].nodes[k]
        pdofe[k] = nodes[n].pdof
        pdofi = nodes[n].pdof
        if pdofi > 0
          Pe[k] = U[pdofi]
        else
          Pe[k] = nodes[n].p
        end
      end
      rho = material[mat][1]
      T6_KeNS!(Ke, De, FeU, FeP, Fe0, Xe, Ue, Pe, rho, nu, coeffnorm)
      for i = 1:12
        dofi = udofe[i]
        if dofi > 0
          F[dofi] -= FeU[i]
          for j = 1:12
            dofj = udofe[j]
            if dofj > 0
              K[dofi, dofj] += Ke[i, j]
            end
          end
          for j = 1:3
            dofj = pdofe[j]
            K[dofi, dofj] += De[i, j]
          end
        end
      end

      for i = 1:3
        dofi = pdofe[i]
        F[dofi] -= FeP[i]
        for j = 1:12
          dofj = udofe[j]
          if dofj > 0
            K[dofi, dofj] += De[j, i]
          end
        end
      end
    end
  end
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
