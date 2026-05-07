function assembly_quad_NS!(P::Parametrisation, pos::Int64, W1::Vector{ComplexF64}, W2::Vector{ComplexF64},
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
            U1e[k, idof] = W1[dof]
            U2e[k, idof] = W2[dof]
          else
            U1e[k, idof] = 0.0
            U2e[k, idof] = 0.0
          end
        end
      end
      nu1 = W1[info.neq+1]
      nu2 = W2[info.neq+1]

      T6_quad_NS!(Fe, Xe, U1e, U2e, nu1, nu2, coeffnorm)

      for i in 1:12
        dofU = dofe[i]
        if dofU > 0
          P.R[dofU, pos] -= Fe[i]
        end
      end
    end

  end  # loop over NE
end

function T6_quad_NS!(Fe::Vector{ComplexF64}, Xe::Matrix{Float64},
  U1e::Matrix{ComplexF64}, U2e::Matrix{ComplexF64}, nu1::ComplexF64, nu2::ComplexF64, coeffnorm::Float64)

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

    Fe[:] += B * (nu1 * [GU2[1, 1]; GU2[2, 2]; 0.5 * (GU2[1, 2] + GU2[2, 1])] + nu2 * [GU1[1, 1]; GU1[2, 2]; 0.5 * (GU1[1, 2] + GU1[2, 1])]) * coeff / coeffnorm
    Fe[:] += 0.5 * N * ([GU1[1, 1] * U2[1] + GU1[1, 2] * U2[2]; GU1[2, 1] * U2[1] + GU1[2, 2] * U2[2]] +
                        [GU2[1, 1] * U1[1] + GU2[1, 2] * U1[2]; GU2[2, 1] * U1[1] + GU2[2, 2] * U1[2]]) * coeff

  end  # loop gauss points

  return nothing

end