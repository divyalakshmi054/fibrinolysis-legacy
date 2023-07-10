# load the include -
include("Include.jl")

#=========== set Visit number in line 20; set [tPA] in line 45 ============#

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis_Interpolation.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Fibrinolysis-Legacy-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 1;
visit_df = filter(:visitid => x->(x==visit), training_df) 

#= # thrombin data 
thrombin_df = CSV.read(joinpath(pwd(),"data","SIM-Actual-P1.csv"),DataFrame)
FIIa = thrombin_df[:,:FIIa]
global FIIa_itp = interpolate(FIIa,BSpline(Linear())) =#

# prothrombin -
FII = visit_df[:,:II]

# size of training set -
(R,C) = size(visit_df)

# time -
t_df = CSV.read(joinpath(pwd(),"actual_legacy_thrombin","SIM-Actual-visit-$(visit)-P1.csv"),DataFrame)
Time = t_df[:,:T]

global FII = visit_df[:,:II]

# main simulation -
SF = 1e9
for i ∈ 1:1

    # build new model -
    dd = deepcopy(model)

    # thrombin data -
    thrombin_df = CSV.read(joinpath(pwd(),"actual_legacy_thrombin","SIM-Actual-visit-$(visit)-P$(i).csv"),DataFrame)
    global FIIa = thrombin_df[:,:FIIa]
    global FIIa_itp = interpolate(FIIa,BSpline(Linear()))

    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = 4.0                    # 1 tPA   SET tPA conc here!
    sfa[2] = visit_df[i,:PAI1]       # 2 PAI1; calculated from literature
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
    α[1] = 0.5
    α[2] = 0.1
    α[3] = 0.015

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
    G[FIIa_idx, 1] = 0.5
    G[FI_idx, 1] = 1.75

    # adjusting parameters for r2
    G[tPA_idx,2] = 0.75    
    # G[Plgn_idx,2] = 0.9
    G[PAI1_idx,2] = -1.0

    # adjusting parameters for r3
    G[Plasmin_idx,3] = 0.75    
    G[TAFI_idx,3] = -1 
    G[FIa_idx,3] = 0.9

    # run the model -
    global (T,U) = evaluate_w_delay(dd,tspan=(0.0,180.0))
    data = [T U]
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,2],sfa[1],FIIa,FII[i])
    print("AUC:" , integrate(T,CF),"\n")
    print("Max: ",maximum(CF),"\n")
    idx_MCF = findfirst(x->x==maximum(CF),CF)
    MCF = maximum(CF)
    if(MCF>=20.0)
        idx_CFT = findfirst(x->x>=20.0,CF)
    elseif(MCF<20.0)
        idx_CFT = findfirst(x->x==MCF,CF)
    end
    idx_2 = findfirst(x->x>=2.0,CF)
    slope = (CF[idx_CFT]-2.0)/(T[idx_CFT] - T[idx_2])
    alpha = atan(slope)*(180/pi)
    print("alpha: ",alpha,"\n")
    
    # dump -
    _PATH_TO_TMP = joinpath(pwd(),"tmp")
    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-visit-$(visit)-Fib-$(tpa_int)-nM-tPA-run-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data,CF),header=vcat("Time",dd.list_of_dynamic_species,"CF")))

    # figures -
    _PATH_TO_FIGS = joinpath(pwd(),"figs")
    path_to_CFfigs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_CF_run$(i).png")
    path_to_thrombin_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_thrombin_run$(i).png")
    #path_to_fibrin_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_fibrin_run$(i).png")
    # path_to_CF_ensemble_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_CF_runs.png")

    Plots.savefig(Plots.plot(T, CF, xticks=0.0:10:180, xlabel="Time (min)", ylabel="CF (mm)", title="Clot firmness vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_CFfigs)
    # Plots.savefig(Plots.plot(T, U[:,3], xticks=0.0:10:180,xlabel="Time (min)", ylabel="FIIa (nM)", title="[Thrombin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_thrombin_figs)
    # Plots.savefig(Plots.plot(T, U[:,4], xticks=0.0:10:180, xlabel="Time (min)", ylabel="FIa (nM)", title="[Fibrin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_fibrin_figs)
    
end