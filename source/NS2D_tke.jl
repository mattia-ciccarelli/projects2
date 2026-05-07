using MAT
using LinearAlgebra
using Statistics
using Plots

function TKE()
    
    param_file = matopen("/home/mattia.ciccarelli/Documents/projects/output/param.mat")
    mappings_r = read(param_file, "mappings_r")
    Avector = read(param_file, "Avector")
    M = read(param_file, "M")
    close(param_file)

    uchannel_file = matopen("/home/mattia.ciccarelli/Documents/results/Uchannel1.mat")
    xlcc = read(uchannel_file, "xlcc")
    close(uchannel_file)

    # 2. Parameters
    LU = 1.0
    coeffnorm = 100.0
    nu0 = coeffnorm * LU / 5200.0 #???????????????? CHANGE
    order_trunc = 7 

    nPages, nNodes, nComponents = size(mappings_r)
    nred = size(Avector, 2) - 1

    # dropdims instead of squeeze 
    condition1 = sum(Avector[:, 1:nred+1], dims=2) .<= order_trunc
    term2use = dropdims(condition1, dims=2) 

    po_all = xlcc
    npo = size(po_all, 2)
    lpo = size(po_all, 1)
    lpo = div(lpo - 2, nred)

   
    Re_list = zeros(Float64, npo)
    TKE_mean = zeros(Float64, npo)

    M_scalar = M[1:nNodes, 1:nNodes]

    println("TKE evaluation with threads()) thread)...")

    # Parallelized cycle on Re number
    Threads.@threads for i in 1:npo
        
        
        nu = po_all[end, i]
        Re_list[i] = coeffnorm * LU / (nu0 + nu)
        
        vars = zeros(Float64, nred + 1)
        vars[end] = nu
        
        U_fluct_time = zeros(Float64, lpo, nNodes, nComponents)
        
        # for cycle on the period
        for j in 1:lpo
            for k in 1:nred
                vars[k] = po_all[k + (j-1)*nred, i]
            end
            
            #monomials
            zprod_all = prod(vars' .^ Avector, dims=2)
            zprod_all_1D = dropdims(zprod_all, dims=2) .* term2use
            
            # Reshape per allineare le dimensioni al tensore 3D
            zprod_all_3D = reshape(zprod_all_1D, (nPages, 1, 1))
            
            # Somma lungo la prima dimensione (pages)
            U_fluct_time[j, :, :] = dropdims(sum(mappings_r .* zprod_all_3D, dims=1), dims=1)
        end
        
        # To take the flactuation it is needed to substitue the mean flow
        #from the total flow field
        U_mean = dropdims(mean(U_fluct_time, dims=1), dims=1)
        TKE_time = zeros(Float64, lpo)
        
        # Spatial integration
        for j in 1:lpo
            # views to avoid the store inside the RAM
            @views u_turb = U_fluct_time[j, :, :] .- U_mean
            @views u_comp = u_turb[:, 1]
            @views v_comp = u_turb[:, 2]
            
        
            TKE_time[j] = 0.5 * (dot(u_comp, M_scalar, u_comp) + dot(v_comp, M_scalar, v_comp))
        end
        #time average 
        TKE_mean[i] = mean(TKE_time)
    end

    println("Evaluation completed")

    # 4. PLOT DEL DIAGRAMMA DI BIFORCAZIONE
    p = plot(Re_list, TKE_mean, 
             linecolor=:blue, 
             linewidth=2, 
             xlabel="Re", 
             ylabel="< TKE >", 
             title="Turbulent Kinetic Energy (ROM Order $order_trunc)",
             xlims=(4000, 6000),
             legend=false,
             grid=true)
    
    display(p)

    return Re_list, TKE_mean
end

@time Re_list, TKE_mean = TKE();