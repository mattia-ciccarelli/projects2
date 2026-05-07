function dpim(K::SparseMatrixCSC{Float64}, M::SparseMatrixCSC{Float64},
  D::Vector{ComplexF64}, VR::Matrix{ComplexF64}, VL::Matrix{ComplexF64}, VR_nu::Vector{Float64},
  nodes::Vector{Snode}, T6::Vector{ST6}, B3::Vector{SB3})

  info.nmm = length(info.Lmm)   # master modes
  nmm = info.nmm

  neq = info.neq
  info.nza = 2 * nmm
  nA = info.nA

  info.nMat = nA + info.nza  # dim of system to be solved

  println("Init Parametrisation")
  P = initParametrisation!(info)

  println("Init System")
  Rhs = Array{ComplexF64}(undef, info.nMat)
  Sol = Array{ComplexF64}(undef, info.nMat)
  Mat = spzeros(ComplexF64, info.nMat, info.nMat)

  println("Order 1")
  #  omega0 = zeros(Float64,nmm)
  BY = Array{ComplexF64}(undef, nA, info.nza)
  XTB = Array{ComplexF64}(undef, info.nza, nA)

  # first eigenval is 0
  for i = 1:nmm
    mm = info.Lmm[i]
    mmcj = info.Lmmcj[i]
    P[1].f[i, i] = D[mm]
    P[1].f[i+nmm, i+nmm] = D[mmcj]
    P[1].W[1:neq, i] = VR[:, mm]
    P[1].W[1:neq, i+nmm] = VR[:, mmcj]
    BY[1:neq, i] = M * VR[1:neq, mm]
    BY[1:neq, i+nmm] = M * VR[1:neq, mmcj]
    XTB[i, 1:neq] = transpose(VL[1:neq, mm]) * M
    XTB[i+nmm, 1:neq] = transpose(VL[1:neq, mmcj]) * M
  end

  # places autovett=0 in last position
  P[1].W[1:nA, info.nza+1] = VR_nu[:]
  P[1].W[nA+1, info.nza+1] = 1

  Lambda = Vector{ComplexF64}(undef, info.nza + 1)
  for i in 1:info.nza
    Lambda[i] = P[1].f[i, i]
  end
  Lambda[info.nza+1] = 0.0

  # no nonautonomous

  println("Higher orders")
  for p = 2:info.max_order
    println("Order $p")
    fillrhs_quad!(nodes, T6, B3, P, p)
    fillWf!(P, p)
    for i in 1:P[p].m # for every alpha vector
      corresp = P[p].corresp[i]
      if corresp > 0
        #        fillWfnonaut!(P,p,P[p].Av[i],i)
        homological!(Sol, Rhs, Mat, Lambda, P[p], i, K, M, BY, XTB)
        P[p].analysed[i] = 1
      elseif corresp < 0
        P[p].W[:, i] = conj(P[p].W[:, -corresp])
        P[p].f[1:nmm, i] = conj(P[p].f[nmm+1:2*nmm, -corresp])
        P[p].f[nmm+1:2*nmm, i] = conj(P[p].f[1:nmm, -corresp])
        P[p].analysed[i] = 1
      end
    end

  end

  return P

end



function homological!(Sol::Vector{ComplexF64}, Rhs::Vector{ComplexF64}, Mat::SparseMatrixCSC{ComplexF64},
  Lambda::Vector{ComplexF64}, P::Parametrisation, Apos::Int64,
  K::SparseMatrixCSC{Float64}, M::SparseMatrixCSC{Float64},
  BY::Matrix{ComplexF64}, XTB::Matrix{ComplexF64})


  uneq = info.neq
  nA = info.nA
  nza = info.nza

  σ = dot(P.Av[Apos], Lambda[1:nza+1])

  resonant_modes = zeros(Bool, info.nza + 1)
  if info.style == 'c'
    LambdaSym = info.LambdaSym
    σ_sym = dot(P.Av[Apos], LambdaSym)
    fill!(resonant_modes, false)
    for i = 1:info.nza+1
      λ_sym = LambdaSym[i]
      if (σ_sym == λ_sym)
        resonant_modes[i] = true
      end
    end
  elseif info.style == 'g'
    fill!(resonant_modes, true)
  else
    error("style not recognized")
  end
  println(P.Av[Apos], "  ", resonant_modes)


  fill!(Rhs, 0.0)
  Rhs[1:nA] = P.R[:, Apos]
  Rhs[1:nA] -= M * P.Wf[1:nA, Apos]

  fill!(Mat.nzval, 0.0)
  Mat[1:uneq, 1:uneq] = σ * M + K

  for j = 1:info.nza
    if resonant_modes[j]  # if resonant
      Mat[1:nA, nA+j] = BY[:, j]
      Mat[nA+j, 1:nA] = XTB[j, :]
    else
      Mat[nA+j, nA+j] = 1
    end
  end

  ps = MKLPardisoSolver()
  solve!(ps,Sol,Mat,Rhs)
  #Sol = Mat \ Rhs

  P.W[1:nA, Apos] = Sol[1:nA]
  P.f[1:info.nza, Apos] = Sol[nA+1:nA+info.nza]

end


function fillrhs_quad!(nodes::Vector{Snode}, T6::Vector{ST6}, B3::Vector{SB3},
  P::Vector{Parametrisation}, p::Int64)

  for p1 in 1:p-1
    p2 = p - p1
    for k1 in 1:P[p1].m, k2 in 1:P[p2].m
      Av = P[p1].Av[k1] + P[p2].Av[k2]
      pos = findfirst(x -> x == Av, P[p].Av)
      if P[p].corresp[pos] > 0
        W1 = P[p1].W[:, k1] # second part is U+psi, first is V
        W2 = P[p2].W[:, k2]
        assembly_quad_NS!(P[p], pos, W1, W2, nodes, T6)
      end
    end
  end

end

function fillWf!(P::Vector{Parametrisation}, p::Int64)

  for p1 in 2:p-1, p2 in 2:p-1
    if (p1 + p2) == p + 1
      for i in 1:P[p1].m
        Av1 = P[p1].Av[i][:]
        for j in 1:P[p2].m
          Av2 = P[p2].Av[j][:]
          for s in 1:info.nza+1
            if Av1[s] > 0
              Av = Av1 + Av2
              Av[s] -= 1
              pos = findfirst(x -> x == Av, P[p].Av)
              P[p].Wf[:, pos] += Av1[s] * P[p1].W[:, i] * P[p2].f[s, j]
            end
          end
        end
      end
    end
  end

end
