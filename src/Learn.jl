# it's training time -

function learn_optim(index::Int, model::BSTModel, training_df::DataFrame; 
    pₒ::Union{Nothing,Array{Float64,1}} = nothing, FIIa_itp,FII)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 4.0                        # 1 tPA
    sfa[2] = 0.5                        # 2 PAI1; calculated from literature
    sfa[3] = training_df[index,:TAFI]   # 3 TAFI
    sfa[4] = training_df[index,:AT]     # 4 AT
    tPA_value = sfa[1];  
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[index, :Fbgn]    # 2 FI / Fbgn
    xₒ[4] = visit_df[index, :Plgn]    # 4 Plgn
    model.initial_condition_array = xₒ
 

    # what is the output array?
    Y = Array{Float64,1}(undef,4)
    Y[1] = visit_df[index, :CT]
    Y[2] = visit_df[index, :MCF]
    # Y[3] = visit_df[index, :LOT]
   # Y[4] = visit_df[index, :LT]
    Y[3] = visit_df[index,:alpha]
    Y[4] = visit_df[index, :AUC]

    # setup initial parameter values and bounds array -
    κ = [
            
            # default: hand fit set -
            0.5         0.05 2.0      ; # 1
            0.1        0.01 2.0      ; # 2
            0.015       0.01 2.0      ; # 3
            0.5         0.05 2.0      ; # 4
            1.75         0.5 5.0       ; # 5
            0.75         0.5 2.0       ; # 6
            1.0         0.5 2.0      ; # 7
            1.0         0.5 2.0       ; # 8
            0.75        0.5 2.0       ; # 9
            1.0        0.5 2.0      ; # 10
            0.90        0.5 2.0       ; # 11
        ];

    # set default set as the start -
    if (isnothing(pₒ) == true)
        P = length(κ[:,1])
        σ = 0.1 # move up to 10%
        pₒ = κ[:,1].*(1 .+ σ*rand(-1:1,P))
    end

    # setup the objective function -
    inner_optimizer = NelderMead()
    obj_function(p) =  loss_scalar(p, Y, model,sfa[1],FII,index)
    results = optimize(obj_function, κ[:,2], κ[:,3], κ[:,1], Fminbox(inner_optimizer),
        Optim.Options(time_limit = 600, show_trace = true, show_every = 10, iterations=100))
    
    # grab the best parameters -
    p_best = Optim.minimizer(results)
    
    # run the sim w/the best parameters -
    # 1 - 9 : α vector
    model.α = p_best[1:3]
    G = model.G
    FII_idx = findfirst(x->x=="FII",model.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
    AT_idx = findfirst(x->x=="AT",model.total_species_list)
    FI_idx = findfirst(x->x=="FI",model.total_species_list)
    tPA_idx = findfirst(x->x=="tPA",model.total_species_list)
    Plgn_idx = findfirst(x->x=="Plgn",model.total_species_list)
    PAI1_idx = findfirst(x->x=="PAI1",model.total_species_list)
    Plasmin_idx = findfirst(x->x=="Plasmin",model.total_species_list)
    TAFI_idx = findfirst(x->x=="TAFI",model.total_species_list)
    FIa_idx = findfirst(x->x=="FIa",model.total_species_list)

    # adjusting parameters for r1
    G[FIIa_idx, 1] = p_best[4]
    G[FI_idx, 1] = p_best[5]

    # adjusting parameters for r2
    G[tPA_idx,2] = p_best[6]    
    G[Plgn_idx,2] = p_best[7]
    G[PAI1_idx,2] = -1*p_best[8]

    # adjusting parameters for r3
    G[Plasmin_idx,3] = p_best[9]    
    G[TAFI_idx,3] = -1*p_best[10]  
    G[FIa_idx,3] = p_best[11]

    # run the model -
    (T,U) = evaluate_w_delay(model, tspan = (0.0, 180.0))
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII)
    data = [T U]
    Xₘ = hcat(data,CF)
    Yₘ = model_output_vector(T, CF) # properties of the CF curve 
    
    return (p_best, T, Xₘ, Yₘ, Y)
end