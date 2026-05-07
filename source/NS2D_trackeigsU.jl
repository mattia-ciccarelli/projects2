using Hungarian

function track_eigenvalues_Uchannel(Rei, Ref, n)

  fname = "./input/" * input_file
  include(fname)
  info.mesh_file = file

  T6, B3, nodes = readgmsh!(info)
  K, M, F0 = analysis(nodes, T6, info)

  function custom_sort(eigenvalues::Vector{Complex{T}}) where T
    sorted_indices = sortperm(abs.(eigenvalues), rev = true)
    sorted_indices = sortperm(imag.(eigenvalues[sorted_indices]))
    return eigenvalues[sorted_indices]
  end

  Re = Rei
  dRe = (Ref- Rei)/n
  nu0 = material[1][2]
  eigens = []
  iter = 1
  while Re <= Ref
    println("Calculation number: ", iter)
    material[1][2] = coeffnorm / Re
    #material[1][2] = nu0 * Rei/Re
    K, M, F0 = analysis_eigs(nodes, T6, info)
    D, VR, VL, VR_nu = eigAB_Uchannel(K, M, F0, nodes)
    sorted_eigenvalues = custom_sort(D)
    push!(eigens, sorted_eigenvalues)
    Re += dRe
    iter+=1
  end

  max_len = maximum(length.(eigens))
  println("Numero massimo di modi trovati in un singolo step: ", max_len)

  for i in 1:length(eigens)
      diff_len = max_len - length(eigens[i])
      if diff_len > 0
          # Aggiungiamo complex(NaN, NaN) per riempire i buchi
          append!(eigens[i], fill(complex(NaN, NaN), diff_len))
      end
  end

  eigens_in_time = hcat(eigens...)
  N, T = size(eigens_in_time)

  matwrite("./output/reconstructed_series2.mat", Dict("time_series" => eigens_in_time))
  
end

function analysis_eigs(nodes, T6, info)


  F = zeros(Float64, info.neq)
  F0 = zeros(Float64, info.neq)

  for iload in 1:10
    residuum = 100
    iter = 0
    while residuum > 1e-5
      #while iter<5

      iter += 1

      global K   # to make it visible out of the while loop

      # faster and safer to reallocate at every iter, as the number of nozero coeffs
      # may vary from iter to iter
      K = SparseMatrixLNK(info.neq, info.neq)
      if iter > 0
        fill!(F, 0.0)
        fill!(F0, 0.0)
      end

      println("Assembling")

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
          #    T6_KeNL!(Ke,Fe,Xe,Ue,material[mat])
          rho, nu = material[mat][1:2]
          nu = nu*10/iload
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

      # Solution phase
      println("Solving system")
      K = SparseMatrixCSC(K)  # convert to other sparse format
      sol = K \ F
      #ps = MKLPardisoSolver()
      #solve!(ps,sol,K,F)

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
      println("Residuum: ", residuum)

    end # while residuum
  end

  ##############################
  # mass assembler
  ##############################

  M = SparseMatrixLNK(info.neq, info.neq)

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

  return K, M, F0
end

function isapproxcomplex(a::Complex, b::Complex; atol::Real = 1e-8, rtol::Real = 1e-8)
  return isapprox(abs(a), abs(b), atol = atol, rtol = rtol) && isapprox(imag(a), imag(b), atol = atol, rtol = rtol)
end

function complexlessthan(a::Complex, b::Complex)
  abs_a = abs(a)
  abs_b = abs(b)

  # First compare by absolute value
  if abs_a != abs_b
    return abs_a < abs_b
  end

  # Then compare by imaginary part
  imag_a = imag(a)
  imag_b = imag(b)
  if imag_a != imag_b
    return imag_a < imag_b
  end

  # Finally compare by real part
  return real(a) < real(b)
end
