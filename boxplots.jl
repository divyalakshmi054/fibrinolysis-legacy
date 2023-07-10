# include -
include("Include.jl")
using StatsPlots

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

# let us plot -
sim_df = DataFrame([MCF_sim],[:simulated])
exp_df = DataFrame([MCF_exp],[:actual])
df = hcat(exp_df,sim_df)
stacked_df = stack(df,[:actual,:simulated])
@df stacked_df boxplot(string.(:variable), :value, fillalpha=0.75, linewidth=2,label="",color=colorant"#6A6A6A")
y = ylabel!("Maximum Clot Firmness (mm)")
savefig(y,joinpath(pwd(),"figs","MCF_V$(visit).pdf"))
savefig(y,joinpath(pwd(),"figs","MCF_V$(visit).png"))