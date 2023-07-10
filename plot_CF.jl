# include -
include("Include.jl")
using Statistics
size_t = 18001;
_PATH_TO_ACTUAL_ENSEMBLE = joinpath(pwd(),"ensemble_redo")

visit_1_CF = Array{Float64}(undef,(size_t,1))
sd_visit_1_CF = Array{Float64}(undef,(size_t,1))
visit_6_CF = Array{Float64}(undef,(size_t,1))
sd_visit_6_CF = Array{Float64}(undef,(size_t,1))

SIM_CF = Array{Float64}(undef,(size_t,10))
Z = zeros(size_t)

T_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Actual-P1-visit-1-tPA-4-nM.csv"),DataFrame)
Time = T_df[1:size_t,:Time]
for i ∈ 1:10
    CF_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Actual-P$(i)-visit-1-tPA-4-nM.csv"),DataFrame)
    Temp_CF = CF_df[1:size_t,:CF]
    SIM_CF[:,i] = Temp_CF
end

visit_1_CF = mean(SIM_CF,dims = 2)
sd_visit_1_CF = Statistics.std(SIM_CF,dims = 2)
U = sd_visit_1_CF
plot(Time,visit_1_CF,color=colorant"#0077BB",label="",lw=2,fillrange=(visit_1_CF-U,U+visit_1_CF),fillalpha=0.4)
xlabel!("Time (min)")
ylabel!("Clot Firmness (mm)")
plot!([-5],[-5],xlim=(0.0,180.0),ylim=(0.0,27.0),line=:scatter,markerstrokecolor=colorant"#0077BB",color=colorant"#0077BB",label = "Visit 1",foreground_color_legend = nothing)

for i ∈ 1:10
    CF_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"SIM-Actual-P$(i)-visit-3-tPA-4-nM.csv"),DataFrame)
    Temp_CF = CF_df[1:size_t,:CF]
    SIM_CF[:,i] = Temp_CF
end

visit_6_CF = mean(SIM_CF,dims = 2)
sd_visit_6_CF = std(SIM_CF,dims = 2)
U = sd_visit_6_CF
plot!(Time,visit_6_CF,color=colorant"#EE7733",label="",lw=2,fillrange=(visit_6_CF-U,U+visit_6_CF),fillalpha=0.4)
y = plot!([-5],[-5],xlim=(0.0,180.0),ylim=(0.0,27.0),line=:scatter,color=colorant"#EE7733",markerstrokecolor = colorant"#EE7733",label = "Visit 6",foreground_color_legend = nothing)
savefig(y,joinpath(pwd(),"figs","Actual-N10-CF-New-New.pdf"))