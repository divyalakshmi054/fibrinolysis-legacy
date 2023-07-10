# load the include -
include("Include.jl")

# parameters -
_PATH_TO_ACTUAL_ENSEMBLE = joinpath(pwd(),"actual_ensemble")

# number of patients -
number_of_patients = 12;

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis_Interpolation.bst")

# build the default model structure -
model = build(path_to_model_file)
_PATH_TO_DATA = joinpath(pwd(),"data")

# load the training data -
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
full_df = CSV.read(path_to_training_data, DataFrame)

# let's filter visit 4 -
visit = 4;
tPA = 4;        # for later
visit_df = filter(:Visit => x->(x==visit), full_df)

# prothrombin -
Prothrombin = visit_df[:,:II]

SF = 1e9

for i ∈ 1:number_of_patients

    if (i==5)   # patient 5, visit 4 data was missing
        continue
    end

    pset_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"PSET-Actual-P$(i)-visit-4-tPA-4-nM.csv"),DataFrame)
    # build new model -
    dd = deepcopy(model)

    thrombin_df = CSV.read(joinpath(pwd(),"data","thrombin","SIM-visit-$(visit)-TF-$(i).csv"),DataFrame)
    FIIa = thrombin_df[:,:FIIa]
    global FIIa_itp = interpolate(FIIa,BSpline(Linear()))
    global FII = visit_df[:,:II]

    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = 4.0                    # 1 tPA   SET tPA conc here!
    sfa[2] = 0.5                    # 2 PAI1; calculated from literature
    sfa[3] = visit_df[i,:TAFI]      # 3 TAFI
    sfa[4] = visit_df[i,:AT]        # 4 AT
    sfa[5] = (1e-14)*SF             # 5 FIIa  
    tpa_int = Int64(sfa[1])
    
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = dd.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :Fbgn]    # 2 FI / Fbgn
    xₒ[4] = visit_df[i, :Plgn]    # 4 Plgn
    dd.initial_condition_array = xₒ

    #update α -
    α = dd.α
    α[1] = pset_df[2,:parameters]
    α[2] = pset_df[3,:parameters]
    α[3] = pset_df[4,:parameters]

    #update G -
    G = dd.G
    
    FII_idx = findfirst(x->x=="FII",dd.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",dd.total_species_list)
    AT_idx = findfirst(x->x=="AT",dd.total_species_list)
    FI_idx = findfirst(x->x=="FI",dd.total_species_list)
    tPA_idx = findfirst(x->x=="tPA",dd.total_species_list)
    Plgn_idx = findfirst(x->x=="Plgn",dd.total_species_list)
    PAI1_idx = findfirst(x->x=="PAI1",dd.total_species_list)
    Plasmin_idx = findfirst(x->x=="Plasmin",dd.total_species_list)
    TAFI_idx = findfirst(x->x=="TAFI",dd.total_species_list)
    FIa_idx = findfirst(x->x=="FIa",dd.total_species_list)

    # adjusting parameters for r1
    G[FIIa_idx, 1] = pset_df[5,:parameters]
    G[FI_idx, 1] = pset_df[6,:parameters]

    # adjusting parameters for r2
    G[tPA_idx,2] = pset_df[7,:parameters]    
    G[Plgn_idx,2] = pset_df[8,:parameters]
    G[PAI1_idx,2] = pset_df[9,:parameters]

    # adjusting parameters for r3
    G[Plasmin_idx,3] = pset_df[10,:parameters]   
    G[TAFI_idx,3] = -1*pset_df[11,:parameters]    
    G[FIa_idx,3] = pset_df[12,:parameters]

    # run the model -
    global (T,U) = evaluate_w_delay(dd,tspan=(0.0,180.0))
    data = [T U]
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII[i])
    
    # dump -
    _PATH_TO_TMP_VALIDATION = joinpath(pwd(),"tmp_validation")
    path_to_sim_data = joinpath(_PATH_TO_TMP_VALIDATION, "SIM-visit-$(visit)-Fib-$(tpa_int)-nM-tPA-run-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data,CF),header=vcat("Time",dd.list_of_dynamic_species,"CF")))

    # figures -
    _PATH_TO_FIGS_VALIDATION = joinpath(pwd(),"figs_validation")
    path_to_CFfigs = joinpath(_PATH_TO_FIGS_VALIDATION, "tPA_$(tpa_int)nM_visit_$(visit)_CF_run$(i).png")
    path_to_thrombin_figs = joinpath(_PATH_TO_FIGS_VALIDATION, "tPA_$(tpa_int)nM_visit_$(visit)_thrombin_run$(i).png")
    path_to_fibrin_figs = joinpath(_PATH_TO_FIGS_VALIDATION, "tPA_$(tpa_int)nM_visit_$(visit)_fibrin_run$(i).png")
    # path_to_CF_ensemble_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_CF_runs.png")

    Plots.savefig(Plots.plot(T, CF, xticks=0.0:10:180, xlabel="Time (min)", ylabel="CF (mm)", title="Clot firmness vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_CFfigs)
    # Plots.savefig(Plots.plot(T, U[:,3], xticks=0.0:10:180,xlabel="Time (min)", ylabel="FIIa (nM)", title="[Thrombin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_thrombin_figs)
    # Plots.savefig(Plots.plot(T, U[:,4], xticks=0.0:10:180, xlabel="Time (min)", ylabel="FIa (nM)", title="[Fibrin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_fibrin_figs)
    
end
