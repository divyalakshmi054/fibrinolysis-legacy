# load the include -
include("Include.jl")

# build the model structure -
_PATH_TO_MODEL = joinpath(pwd(),"model")
path_to_model_file = joinpath(_PATH_TO_MODEL, "Coagulation.bst")
#model_buffer = read_model_file(path_to_model_file)

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 4

# let's filter visit 4s since we look to train using that visit
visit_df = filter(:Visit => x->(x==visit), training_df) 


# size of training set -
(R,C) = size(visit_df)


# main simulation loop -
SF = 1e9
for i ∈ 1:R

    # build new model -
    dd = deepcopy(model)

    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = training_df[i,:TFPI]       # 1 TFPI
    sfa[2] = training_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

    # grab the multiplier from the data -
    ℳ = dd.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[i, :II]         # 1 FII 
    xₒ[2] = training_df[i, :VII]        # 2 FVII 
    xₒ[3] = training_df[i, :V]          # 3 FV
    xₒ[4] = training_df[i, :X]          # 4 FX
    xₒ[5] = training_df[i, :VIII]       # 5 FVIII
    xₒ[6] = training_df[i, :IX]         # 6 FIX
    xₒ[7] = training_df[i, :XI]         # 7 FXI
    xₒ[8] = training_df[i, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF                  # 9 FIIa
    dd.initial_condition_array = xₒ

    # update α -
    α = dd.α
    α[1] = 0.061
    α[9] = 0.70

    # setup -
    G = dd.G
    
    idx = findfirst(x->x=="FVIIa",dd.total_species_list)
    G[idx, 4] = 0.1

    # what is the index of TRAUMA?
    idx = findfirst(x->x=="AT",dd.total_species_list)
    G[idx, 9] = 0.045

    # what is the index of TFPI?
    idx = findfirst(x->x=="TFPI",dd.total_species_list)
    G[idx, 1] = -0.65

    # run the model -
    global (T,U) = evaluate(dd)
    #X = hcat(U...)
    #data = [T transpose(X)]
    data = [T U]

    # dump -
    _PATH_TO_TMP = joinpath(pwd(),"tmp")

    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-patient-$(i).dat")
    CSV.write(path_to_sim_data, Tables.table(data))
end

