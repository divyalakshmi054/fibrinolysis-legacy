function loss_scalar(κ::Array{Float64,1}, Y::Array{Float64,1},  model::BSTModel, tPA::Float64, FII::Float64,index::Int)

    # get the G matrix -
    G = model.G

    # κ map
    # 1 - 9 : α vector
    model.α = κ[1:3]

    #Gs -
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
    G[FIIa_idx, 1] = κ[4]
    G[FI_idx, 1] = κ[5]

    # adjusting parameters for r2
    G[tPA_idx,2] = κ[6]    
    G[Plgn_idx,2] = κ[7]
    G[PAI1_idx,2] = κ[8]

    # adjusting parameters for r3
    G[Plasmin_idx,3] = κ[9]    
    G[TAFI_idx,3] = -1*κ[10]  
    G[FIa_idx,3] = κ[11]

    # run the model -
    (T,U) = evaluate_w_delay(model, tspan = (0.0, 180.0))
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],tPA,FIIa,FII)
    Yₘ = model_output_vector(T, CF) # properties of the CF curve
    #if (SciMLBase.successful_retcode(sol))
    ϵ = norm((Y .- Yₘ).^2)
    #else 
     #   ϵ = Inf
    #end

    # @info ϵ

    # return -
    return ϵ
end

function loss(κ::Array{Float64,1}, Y::Array{Float64,1},  model::BSTModel)

    # get the G matrix -
    G = model.G

    # κ map
    # 1 - 9 : α vector
    model.α = κ[1:3]

    #Gs -
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
    G[FIIa_idx, 1] = κ[4]
    G[FI_idx, 1] = κ[5]

    # adjusting parameters for r2
    G[tPA_idx,2] = κ[6]    
    G[Plgn_idx,2] = κ[7]
    G[PAI1_idx,2] = -1*κ[8]

    # adjusting parameters for r3
    G[Plasmin_idx,3] = κ[9]    
    G[TAFI_idx,3] = -1*κ[10]  
    G[FIa_idx,3] = κ[11]

    # run the model -
    (T,U) = evaluate_w_delay(model, tspan = (0.0, 180.0))
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII[index])
    data = [T U]
    Xₘ = hcat(data,CF)
    Yₘ = model_output_vector(T, CF) # properties of the CF curve 
    ϵ = norm((Y .- Yₘ).^2)

    # return -
    return ϵ
end