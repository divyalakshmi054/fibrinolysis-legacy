# include -
include("Include.jl")
using StatsPlots
using CategoricalArrays

# what fingerprints do we have?
number_of_samples = 10;
AUC_sim = zeros(number_of_samples)
MCF_sim = zeros(number_of_samples)
CT_sim = zeros(number_of_samples)
alpha_sim = zeros(number_of_samples)

# for the actual exp data -
AUC_exp = zeros(number_of_samples)
MCF_exp = zeros(number_of_samples)
CT_exp = zeros(number_of_samples)
alpha_exp = zeros(number_of_samples)

# finally, pick a visit -
visit = 3;

# open up our files -
for i ∈ 1:number_of_samples
    df = CSV.read(joinpath(pwd(),"ensemble_redo","OUT-Actual-P$(i)-visit-$(visit)-tPA-4-nM.csv"),DataFrame)
    CT_exp[i] = df[1,:actual]
    CT_sim[i] = df[1,:simulated]
    MCF_exp[i] = df[2,:actual]
    MCF_sim[i] = df[2,:simulated]
    alpha_exp[i] = df[3,:actual]
    alpha_sim[i] = df[3,:simulated]
    AUC_exp[i] = df[4,:actual]
    AUC_sim[i] = df[4,:simulated]
end

nam = (repeat(1:10, outer=2))
#nam = repeat(["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"], outer =2)
sx = repeat(["Simulated","Actual"],inner=10)
y = vcat(alpha_sim,alpha_exp)
toplot = groupedbar(nam,y,group=sx,ylabel="Alpha Angle",label = "",bar_width = 0.67, 
lw = 0, c = [colorant"#1b1b1b" colorant"#a2a2a2"], markerstrokewidth = 1.5, grid = false,xticks=1:1:10,xlim = (0,11), ylim = (0,85))
xlabel!("Patient Index")
scatter!([-5,-5],color=colorant"#1b1b1b",markerstrokecolor=colorant"#1b1b1b",label = "Simulated", foreground_color_legend = nothing, legend = :topright)
scatter!([-5,-5],color=colorant"#a2a2a2",markerstrokecolor=colorant"#a2a2a2",label = "Actual", foreground_color_legend = nothing, legend = :topright)
savefig(toplot,joinpath(pwd(),"figs","Alpha_V$(visit)_BarPlot.pdf"))
savefig(toplot,joinpath(pwd(),"figs","Alpha_V$(visit)_BarPlot.png"))