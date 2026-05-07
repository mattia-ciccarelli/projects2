
function analysis(nodes::Vector{Snode}, T6::Vector{ST6}, info::Sinfo, verbose::Bool=true)

  if !verbose
    println("Solving static solution")
  end
  # assignes inflow
  for n in 1:info.NN
    for d in 1:2
      dof = nodes[n].udof[d]
      if dof == -2
        nodes[n].u[d] = inflow(nodes[n].coor, d)
      end
    end
  end

  if first_run
    uneq = 0
  end
  # equation numbering: disp
  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:6
        n = T6[e].nodes[k]
        for d in 1:2
          if nodes[n].udof[d] == 0
            uneq += 1
            nodes[n].udof[d] = uneq
          end
        end
      end
    end
  end

  if first_run
    pneq = 0
  end
  # equation numbering: press
  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:3
        n = T6[e].nodes[k]
        if nodes[n].pdof == 0
          pneq += 1
          nodes[n].pdof = uneq + pneq
        end
      end
    end
  end
  if first_run
    info.uneq = uneq
    info.pneq = pneq
    info.nK = uneq
    info.nA = uneq + pneq
    neq = uneq + pneq
    info.neq = uneq + pneq
  end
  neq = info.neq
  uneq = info.uneq
  pneq = info.pneq
  sol = zeros(Float64, neq)
  F = zeros(Float64, neq)
  F0 = zeros(Float64, neq)

  for iload in 1:10
    residuum = 100
    iter = 0
    while residuum > 1e-8
      #while iter<5

      iter += 1

      global K   # to make it visible out of the while loop

      # faster and safer to reallocate at every iter, as the number of nozero coeffs
      # may vary from iter to iter
      K = SparseMatrixLNK(neq, neq)
      if iter > 0
        fill!(F, 0.0)
        fill!(F0, 0.0)
      end

      if verbose
      println("Assembling")
      end

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
            Ue[k, :] = nodes[n].u[:]
          end
          for k in 1:3
            n = T6[e].nodes[k]
            pdofe[k] = nodes[n].pdof
            Pe[k] = nodes[n].p
          end
          rho, nu = material[mat][1:2]
          nu = nu/iload*10
          T6_KeNS!(Ke, De, FeU, FeP, Fe0, Xe, Ue, Pe, rho, nu, coeffnorm)

          for i = 1:12
            dofi = udofe[i]
            if dofi > 0
              F[dofi] -= FeU[i]
              F0[dofi] -= Fe0[i]
              for j = 1:12
                dofj = udofe[j]
                if dofj > 0
                  K[dofi, dofj] += Ke[i, j]
                end
              end
              for j = 1:3
                dofj = pdofe[j]
                if dofj > 0
                  K[dofi, dofj] += De[i, j]
                end
              end
            end
          end

          for i = 1:3
            dofi = pdofe[i]
            if dofi > 0
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

      # Solution phase
      if verbose
      println("Solving system")
      end
      K = SparseMatrixCSC(K)  # convert to other sparse format
      K = K + 1e-12 * I
      #sol = K \ F
      ps = MKLPardisoSolver()
      solve!(ps,sol,K,F)

      # fill dofs
      for n in 1:info.NN
        for d in 1:2
          dof = nodes[n].udof[d]
          if dof > 0
            nodes[n].u[d] += sol[dof]
          end
          dof = nodes[n].udof[d]
        end
        dof = nodes[n].pdof
        if dof > 0
          nodes[n].p += sol[dof]
        end
      end

      residuum = norm(sol)
      if verbose
      println("Residuum: ", residuum)
      end

    end # while residuum
  end
  ##############################
  # mass assembler
  ##############################

  M = SparseMatrixLNK(neq, neq)

  # allocation of elements arrays
  Xe = zeros(Float64, 6, 2)
  dofe = zeros(Int64, 12, 1)
  Me = zeros(Float64, (12, 12))

  for e in 1:info.NE
    mat = T6[e].mat
    if mat > 0
      for k in 1:6
        n = T6[e].nodes[k]
        Xe[k, :] = nodes[n].coor
        dofe[2*k-1:k*2] = nodes[n].udof
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

  ##############################
  # GMSH output of static solution
  ##############################

  outgmsh(info.NN, nodes, T6, "_u0")

  return K, M, F0

end


