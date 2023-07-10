# load include -
include("Include.jl")

# open real data -
visit = 4;
tPA = 0;

_PATH_TO_REAL_DATA = joinpath(pwd(),"data","sorted")
real_df = CSV.read(joinpath(_PATH_TO_REAL_DATA,"REAL-visit-$(visit)-TEG-$(tPA)-nM-tPA.csv"),DataFrame)
Time = real_df[:,:Time_sec]

_PATH_TO_SIM_DATA = joinpath(pwd(),"tmp_validation")

number_of_patients = 12;

RMSE = Array{Float64}(undef,(12,1))
# fitness = Array{Float64}(undef,(12,1))

for i ∈ 1:number_of_patients
    if(i==5)
        continue
    end
    sim_df = CSV.read(joinpath(_PATH_TO_SIM_DATA,"SIM-visit-$(visit)-Fib-$(tPA)-nM-tPA-run-$(i).csv"),DataFrame)
    CF_data = sim_df[:,:CF]
    CF_itp = interpolate(CF_data,BSpline(Linear()))
    CF_use = Array{Float64}(undef,(length(Time),1))
    for j ∈ 1:length(Time)
        CF_use[j] = CF_itp(Time[j]*100+1)
    end
    # fitness[i] = norm((CF_use.-real_df[:,2*i]).^2)
    RMSE[i] = sqrt(mse(CF_use,real_df[:,2*i]))
end


