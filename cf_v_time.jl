# include -
include("Include.jl")

# SIM-visit-$(visit)-Fib-$(tpa_int)-nM-tPA-run-$(i).csv

visit = 2;
tPA_int = 4;
patient_ID = ["Patient 303 CF(mm)","Patient 304 CF(mm)","Patient 307 CF(mm)","Patient 309 CF(mm)","Patient 313 CF(mm)","Patient 314 CF(mm)","Patient 316 CF(mm)","Patient 317 CF(mm)","Patient 319 CF(mm)","Patient 320 CF(mm)","Patient 321 CF(mm)","Patient 323 CF(mm)"]

sim_data = CSV.read(joinpath(pwd(),"tmp","SIM-visit-$(visit)-Fib-$(tPA_int)-nM-tPA-run-1.csv"),DataFrame)
time = sim_data[:,:Time]
CF_combined = Array{Float64}(undef,(length(time),12))

for i ∈ 1:12
    sim_data = CSV.read(joinpath(pwd(),"tmp","SIM-visit-$(visit)-Fib-$(tPA_int)-nM-tPA-run-$(i).csv"),DataFrame)
    CF_combined[:,i] = sim_data[:,:CF]
end

path_to_cf_data = joinpath(pwd(),"CF_data","CF-v-time-visit-$(visit)-Fib-$(tPA_int)-nM-tPA.csv")
CSV.write(path_to_cf_data, Tables.table(hcat(time,CF_combined)),header=vcat("Time",patient_ID))