# load the include -
include("Include.jl")

# build performance function -
function performance(κ, model::BSTModel, visit_df::DataFrame, i::Int64, IIa, II)

    # main simulation -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 4.0                    # 1 tPA   SET tPA conc here!
    sfa[2] = 0.5                    # 2 PAI1; calculated from literature
    sfa[3] = visit_df[i,:TAFI]      # 3 TAFI
    sfa[4] = visit_df[i,:AT]        # 4 AT
    sfa[5] = (1e-14)*SF             # 5 FIIa
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :Fbgn]    # 2 FI / Fbgn
    xₒ[4] = visit_df[i, :Plgn]    # 4 Plgn
    model.initial_condition_array = xₒ
    
    #get the parameters -
    tmp_alpha = κ[1:3]
    g = κ[4:end]

    # set new parameters -
    α = model.α
    α = tmp_alpha
    model.α = α

    # set G values -
    G = model.G;

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
    G[FIIa_idx, 1] = g[1]
    G[FI_idx, 1] = g[2]
 
    # adjusting parameters for r2
    G[tPA_idx,2] = g[3]    
    G[Plgn_idx,2] = g[4]
    G[PAI1_idx,2] = -g[5]
 
    # adjusting parameters for r3
    G[Plasmin_idx,3] = g[6]   
    G[TAFI_idx,3] = -g[7]  
    G[FIa_idx,3] = g[8]

    # put it back -
    model.G = G;

    # solve -
    # run the model -
    global (T,U) = evaluate_w_delay(model,tspan=(0.0,180.0))
    data = [T U]
    CF = Array{Float64,1}
    FIIa = IIa;
    FII = II;
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII)

    # test -
    return integrate(T,CF)    # AUC
end

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis_Interpolation.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Fibrinolysis-Legacy-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 3

# let's filter visit 4s since we look to train using that visit
visit_df = filter(:visitid => x->(x==visit), training_df) 

# prothrombin -
II = visit_df[:,:II]

# size of training set -
(R,C) = size(visit_df)

a = [0.5, 0.1, 0.015]

#update G -
# G = model.G
g = [0.5, 1.75, 0.75, 1.0, 1.0, 0.75, 1.0, 0.9] # look at sample_ensemble.jl for specific G values

# fusion -
parameters = vcat(a,g)
low = minimum(parameters)
high = maximum(parameters)

#= norm_params = Array{Float64}(undef,(11,1))

for i ∈ 1:length(parameters)
    norm_params[i] = parameters[i]/high
end

parameters = norm_params =#

np = length(parameters)

L = zeros(np)
U = zeros(np)
for pᵢ ∈ 1:(np)
    L[pᵢ] = 0.5*parameters[pᵢ]
    U[pᵢ] = 2.0*parameters[pᵢ]
end
#L[end] = -3.0;
#U[end] = 0.0;

patient_index = 3;
samples = 10000;
bootreps = 100;
prothrombin = II[patient_index]

# thrombin data -
thrombin_df = CSV.read(joinpath(pwd(),"actual_legacy_thrombin","SIM-Actual-visit-$(visit)-P$(patient_index).csv"),DataFrame)
global FIIa = thrombin_df[:,:FIIa]
global FIIa_itp = interpolate(FIIa,BSpline(Linear()))

 # initialize -
 sampler = SobolSample()
    
 # generate a sampler -
 (A,B) = QuasiMonteCarlo.generate_design_matrices(samples,L,U,sampler)

# setup call to Sobol method -
F(parameters) =  performance(parameters, model, visit_df, patient_index, FIIa, prothrombin)
m = gsa(F,Sobol(order = [0,1], nboot = bootreps, conf_level = 0.95), A,B)

# dump -
results_array = hcat(m.ST, m.S1, m.ST_Conf_Int, m.S1_Conf_Int)

# write -
#CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-second-order-boot-$(bootreps).csv"), Tables.table(m.S2))
#CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-second-order-CI-boot-$(bootreps).csv"), Tables.table(m.S2_Conf_Int))
CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-boot-$(bootreps).csv"), Tables.table(results_array), header = vcat("Total_order", "First_order", "Total_order_CI", "First_order_CI"))