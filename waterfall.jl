# load the include -
include("Include.jl")

# how many samples?
number_of_samples = 12
number_of_parameters = 11
ensemble_archive = zeros(number_of_parameters+1,1); # first row is the fitness 

_PATH_TO_DATA = joinpath(pwd(),"data")

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Fibrinolysis-Legacy-Transformed-w-Labels.csv")
full_df = CSV.read(path_to_training_data, DataFrame)

visit = 1;
tPA = 4        # for later
visit_df = filter(:visitid => x->(x==visit), full_df)
(R,C) = size(visit_df)

# prothrombin -
Prothrombin = visit_df[:,:II]

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis_Interpolation.bst")

# build the default model structure -
model = build(path_to_model_file)

# main loop -
p_previous = nothing
for i ∈ 8:10

    FII = Prothrombin[i]

    # thrombin data -
    thrombin_df = CSV.read(joinpath(pwd(),"actual_legacy_thrombin","SIM-Actual-visit-$(visit)-P$(i).csv"),DataFrame)
    global FIIa = thrombin_df[:,:FIIa]
    global FIIa_itp = interpolate(FIIa,BSpline(Linear()))
    
    # run the learn routine -
    (p, T, Xₘ, Yₘ, Y) = learn_optim(i, model, visit_df; pₒ = nothing, FIIa_itp,FII)

    # compute the fitness -
    fitness = norm((Yₘ .- Y).^2)
    global ensemble_archive[1] = fitness
    
    # cache the parameters -
    for k ∈ 1:number_of_parameters
        global ensemble_archive[k+1] = p[k]
    end

    _PATH_TO_ACTUAL_ENSEMBLE = joinpath(pwd(),"ensemble_redo")

    # dump parameters to disk (just in case) -
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "PSET-Actual-P$(i)-visit-$(visit)-tPA-$(tPA)-nM.csv"), 
        Tables.table(ensemble_archive); header=["parameters"])
    
    # dump output to disk -
    data_output = [Y Yₘ]
    data_output_header = ["actual", "simulated"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "OUT-Actual-P$(i)-visit-$(visit)-tPA-$(tPA)-nM.csv"), 
        Tables.table(data_output); header = data_output_header)

    # dump state to disk -
    # data_state = transpose(vcat(reshape(T,1,length(T)), Xₘ))
    data_state_header = ["Time","FI","FIa","FDP","Plgn","Plasmin","CF"]
    CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Actual-P$(i)-visit-$(visit)-tPA-$(tPA)-nM.csv"),  
        Tables.table(Xₘ); header=vcat("Time",model.list_of_dynamic_species,"CF"))

    # update/clean up for the next patient -
    global p_previous = p
end