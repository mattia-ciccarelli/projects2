

function eigAB(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},F0::Vector{Float64},nodes::Vector{Snode})

println("Computing eigenvalues")

neq=info.neq
neig=info.neig
uneq=info.uneq

mat"""
[$VR,$DR] = eigs(-$K,$M,$neig,0-2.5i);
$DR=diag($DR);
"""
for i=1:neig
 println(i," ",DR[i])
end

mat"""
[$VL,$DL] = eigs(-transpose($K),transpose($M),$neig,0-2.5i);
$DL=diag($DL);
"""
#println(DL)

# normalisation
for i = 1:neig
  cc = transpose(VL[:,i])*M*VR[:,i]
  # for j = 1:neq
  #   VR[j,i] /= sqrt(cc)
  #   VL[j,i] /= sqrt(cc)
  # end
  VR[:,i] /= norm(VR[:,i], Inf)
  VL[:,i] /= norm(VL[:,i], Inf)
  # outgmsheig(info.NN,nodes,real(VR[:,i]),"_eig"*string(i))  # prints displacements
end

VR_nu=zeros(Float64,neq)
ps = MKLPardisoSolver()
solve!(ps,VR_nu,K,F0)
#VR_nu = K \ F0
# VR_nu /= norm(VR_nu, Inf)
# outgmsheig(info.NN,nodes,VR_nu,"_eig"*string(neig+1))  # prints displacements

return DR,VR,VL,VR_nu

end


function eigABiter(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},F0::Vector{Float64},nodes::Vector{Snode})

  println("Computing eigenvalues")

  neq=info.neq
  neig=info.neig
  uneq=info.uneq

  mat"""
  [$VR,$DR] = eigs(-$K,$M,$neig,'SM');
  $DR=diag($DR);
  """
  for i=1:neig
   println(i," ",DR[i])
  end

  mat"""
  [$VL,$DL] = eigs(-transpose($K),transpose($M),$neig,'SM');
  $DL=diag($DL);
  """
  #println(DL)

  # normalisation
  for i = 1:neig
    cc = transpose(VL[:,i])*M*VR[:,i]  # compl conj
    for j = 1:neq   # arbitrary: scales only VR  (seems to coincide with old approach)
      VR[j,i] /= sqrt(cc)
      VL[j,i] /= sqrt(cc)
    end
  #=   VR[:,i] /= norm(VR[:,i], Inf)
    VL[:,i] /= norm(VL[:,i], Inf) =#
    # outgmsheig(info.NN,nodes,real(VR[:,i]),"_eig"*string(i))  # prints displacements
  end

  VR_nu=zeros(Float64,neq)
  ps = MKLPardisoSolver()
  solve!(ps,VR_nu,K,F0)
  # outgmsheig(info.NN,nodes,VR_nu,"_eig"*string(neig+1))  # prints displacements

  return DR,VR,VL,VR_nu

  end

  function eigAB_Uchannel(K::SparseMatrixCSC{Float64}, M::SparseMatrixCSC{Float64}, F0::Vector{Float64}, nodes::Vector{Snode})

    println("Computing eigenvalues with multiple shifts")

    neq = info.neq
    neig = info.neig
    # uneq = info.uneq 

    shifts = [+2.5im +7.5im -2.5im -7.5im]
    
    DR_all = ComplexF64[]
    VR_all = Matrix{ComplexF64}(undef, neq, 0)
    VL_all = Matrix{ComplexF64}(undef, neq, 0)

    for actuals in shifts
        println("actual shift = ", actuals)

        mat"""
        [$VR_tmp,$DR_tmp] = eigs(-$K,$M,$neig,$actuals);
        $DR_tmp=diag($DR_tmp);
        
        [$VL_tmp,$DL_tmp] = eigs(-transpose($K),transpose($M),$neig,$actuals);
        $DL_tmp=diag($DL_tmp);
        """
        
        append!(DR_all, DR_tmp)
        VR_all = hcat(VR_all, VR_tmp) 
        VL_all = hcat(VL_all, VL_tmp)
    end

    println("Filtering")
    tol = 1e-5 
    indices = Int[]
    
    for i = 1:length(DR_all)
        double = false
        for j in indices
            if abs(DR_all[i] - DR_all[j]) < tol
                double = true
                break
            end
          end
        if !double
            push!(indices, i)
        end
    end

    # Extraction of unique eigenvalues/vector
    DR_unique = DR_all[indices]
    VR_unique = VR_all[:, indices]
    VL_unique = VL_all[:, indices]
    
    num_unique = length(DR_unique)
    
    #normalization
    for i = 1:num_unique
        #cc = transpose(VL_unique[:,i]) * M * VR_unique[:,i]
        
        VR_unique[:,i] /= norm(VR_unique[:,i], Inf)
        VL_unique[:,i] /= norm(VL_unique[:,i], Inf)
    end

    VR_nu = zeros(Float64, neq)
    VR_nu = K \ F0

    println("\nFinal eigenvalue spectrum:")
    for i = 1:num_unique
        println("Mode ", i, ": ", DR_unique[i])
    end

    return DR_unique, VR_unique, VL_unique, VR_nu

end



function eigAB_Uchannel_parall(K::SparseMatrixCSC{Float64}, M::SparseMatrixCSC{Float64}, F0::Vector{Float64}, nodes)

    println("Computing eigenvalues with multiple shifts via Multi-Processing")

    neq = info.neq
    neig = info.neig

    shifts = [+2.5im, +7.5im, -2.5im, -7.5im]
    
    println("Converting matrices to ComplexF64...")
    Kc = ComplexF64.(K)
    Mc = ComplexF64.(M)
    Kc_adj = sparse(Kc')
    Mc_adj = sparse(Mc')

    # ---> MAGIA DISTRIBUITA: pmap gestisce l'invio e la ricezione sicura
    risultati = pmap(shifts) do actuals
        # Questo codice viene eseguito da un Worker completamente isolato
        println("Processo $(myid()) sta calcolando lo shift = $actuals")
        
        DR_tmp, VR_tmp = eigs(-Kc, Mc; nev=neig, sigma=actuals)
        DL_tmp, VL_tmp = eigs(-Kc_adj, Mc_adj; nev=neig, sigma=actuals)
        
        # Il Worker "impacchetta" i risultati in una tupla e li rispedisce al Master
        return (DR_tmp, VR_tmp, VL_tmp)
    end

    # ---> SPACCHETTIAMO I RISULTATI 
    # (risultati è un vettore di tuple ordinato esattamente come gli shifts)
    DR_results = [r[1] for r in risultati]
    VR_results = [r[2] for r in risultati]
    VL_results = [r[3] for r in risultati]

    # 3. UNIONE SICURA (Seriale, sul Master): assembliamo i pezzi
    DR_all = reduce(vcat, DR_results)
    VR_all = reduce(hcat, VR_results)
    VL_all = reduce(hcat, VL_results)

    println("Filtering")
    tol = 1e-5 
    indices = Int[]
    
    # Il resto del tuo codice per filtrare rimane identico!
    for i = 1:length(DR_all)
        double = false
        for j in indices
            if abs(DR_all[i] - DR_all[j]) < tol
                double = true
                break
            end
        end
        if !double
            push!(indices, i)
        end
    end

    # Extraction of unique eigenvalues/vector
    DR_unique = DR_all[indices]
    VR_unique = VR_all[:, indices]
    VL_unique = VL_all[:, indices]
    
    num_unique = length(DR_unique)
    
    # Normalization
    for i = 1:num_unique
        VR_unique[:,i] /= norm(VR_unique[:,i], Inf)
        VL_unique[:,i] /= norm(VL_unique[:,i], Inf)
    end

    # Risoluzione sistema lineare
    VR_nu = K \ F0

    println("\nFinal eigenvalue spectrum:")
    for i = 1:num_unique
        println("Mode ", i, ": ", DR_unique[i])
    end

    return DR_unique, VR_unique, VL_unique, VR_nu
end


function eigAB_Uchannel_parall2(K::SparseMatrixCSC{Float64}, M::SparseMatrixCSC{Float64}, F0::Vector{Float64}, nodes)

    # 1. TRUCCO ANTI-STALLO: Spegniamo i thread di BLAS per evitare conflitti con i thread di Julia
    LinearAlgebra.BLAS.set_num_threads(1)

    println("Computing eigenvalues via KrylovKit (Pure Julia Multi-Threading)")

    neq = info.neq
    neig = info.neig

    shifts = [+2.5im, +7.5im, -2.5im, -7.5im]
    n_shifts = length(shifts)
    
    DR_results = Vector{Vector{ComplexF64}}(undef, n_shifts)
    VR_results = Vector{Matrix{ComplexF64}}(undef, n_shifts)
    VL_results = Vector{Matrix{ComplexF64}}(undef, n_shifts)

    println("Converting matrices to ComplexF64...")
    Kc = ComplexF64.(K)
    Mc = ComplexF64.(M)
    Kc_adj = sparse(Kc')
    Mc_adj = sparse(Mc')

    Threads.@threads for i in 1:n_shifts
        actuals = shifts[i]
        println("Thread ", Threads.threadid(), " computing shift = ", actuals)

        # ---------------------------------------------------------
        # AUTOVETTORI DESTRI
        F_R = -Kc - actuals * Mc
        
        println("Thread ", Threads.threadid(), ": Inizio fattorizzazione LU Destra (può richiedere un po')...")
        fact_R = lu(F_R) 
        println("Thread ", Threads.threadid(), ": Fattorizzazione LU Destra COMPLETATA!")
        
        op_R = v -> fact_R \ (Mc * v)
        v0_R = rand(ComplexF64, neq)

        println("Thread ", Threads.threadid(), ": Inizio iterazioni KrylovKit Destra...")
        # 2. TOLLERANZA E MAXITER: Diciamo a KrylovKit di accontentarsi di 1e-6 e di non fare iterazioni infinite
        vals_inv_R, vecs_R, _ = eigsolve(op_R, v0_R, neig, :LM, krylovdim = neig * 2, tol = 1e-6, maxiter = 300)
        println("Thread ", Threads.threadid(), ": Autovalori Destri TROVATI!")

        # ---------------------------------------------------------
        # AUTOVETTORI SINISTRI (Trasposti)
        F_L = -Kc_adj - actuals * Mc_adj
        
        println("Thread ", Threads.threadid(), ": Inizio fattorizzazione LU Sinistra...")
        fact_L = lu(F_L)
        println("Thread ", Threads.threadid(), ": Fattorizzazione LU Sinistra COMPLETATA!")
        
        op_L = v -> fact_L \ (Mc_adj * v)
        v0_L = rand(ComplexF64, neq)

        println("Thread ", Threads.threadid(), ": Inizio iterazioni KrylovKit Sinistra...")
        # Stessa tolleranza qui
        vals_inv_L, vecs_L, _ = eigsolve(op_L, v0_L, neig, :LM, krylovdim = neig * 2, tol = 1e-6, maxiter = 300)
        println("Thread ", Threads.threadid(), ": Autovalori Sinistri TROVATI!")
        
        # ---------------------------------------------------------
        # 3. IL TAGLIO DI SICUREZZA
        # =========================================================
        # Troviamo la dimensione più piccola tra i destri, i sinistri e il numero richiesto
        n_keep = min(length(vals_inv_R), length(vals_inv_L), neig)
        
        println("Thread ", Threads.threadid(), ": Sto per salvare $n_keep autovalori per il cassetto $i.")

        # Applichiamo il limite [1:n_keep] a TUTTE le estrazioni
        DR_tmp = actuals .+ (1 ./ vals_inv_R[1:n_keep])
        VR_tmp = reduce(hcat, vecs_R[1:n_keep])
        VL_tmp = reduce(hcat, vecs_L[1:n_keep])

        # Salvataggio sicuro e di eguale misura nei cassetti
        DR_results[i] = DR_tmp
        VR_results[i] = VR_tmp
        VL_results[i] = VL_tmp

    # ... [IL RESTO DEL CODICE RIMANE IDENTICO, da DR_all = reduce... in poi] ...
    end

    # ---> ASSEMBLAGGIO SERIALE (Sicuro) <---
    DR_all = reduce(vcat, DR_results)
    VR_all = reduce(hcat, VR_results)
    VL_all = reduce(hcat, VL_results)

    println("Filtering...")
    tol = 1e-5 
    indices = Int[]
    
    # [IL FILTRAGGIO RIMANE IDENTICO AL TUO CODICE ORIGINALE]
    for i = 1:length(DR_all)
        double = false
        for j in indices
            if abs(DR_all[i] - DR_all[j]) < tol
                double = true
                break
            end
        end
        if !double
            push!(indices, i)
        end
    end

    DR_unique = DR_all[indices]
    VR_unique = VR_all[:, indices]
    VL_unique = VL_all[:, indices]
    
    num_unique = length(DR_unique)
    
    for i = 1:num_unique
        VR_unique[:,i] /= norm(VR_unique[:,i], Inf)
        VL_unique[:,i] /= norm(VL_unique[:,i], Inf)
    end

    VR_nu = K \ F0

    println("\nFinal eigenvalue spectrum:")
    for i = 1:num_unique
        println("Mode ", i, ": ", DR_unique[i])
    end

    return DR_unique, VR_unique, VL_unique, VR_nu
end