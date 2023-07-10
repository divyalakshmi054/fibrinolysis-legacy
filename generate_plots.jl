include("Include.jl")

# which visit?
visit = 2

# tpa level?
tpa = 4

# loading up simulation data -
_PATH_TO_TMP_VALIDATION = joinpath(pwd(),"tmp_validation")
# load time vector from a sample simulation file (should all have same times)
time_df = CSV.read(joinpath(_PATH_TO_TMP_VALIDATION,"SIM-visit-$(visit)-Fib-$(tpa)-nM-tPA-run-1.csv"),DataFrame)
t = time_df[:,:Time]

sim_data = Array{Float64}(undef,(length(t),12))   # change number of columns to match number of sims/runs

for i ∈ 1:12
    
    if(i==5)
        continue
    end

    _PATH_TO_FILE = joinpath(_PATH_TO_TMP_VALIDATION,"SIM-visit-$(visit)-Fib-$(tpa)-nM-tPA-run-$(i).csv")
    sim_df = CSV.read(_PATH_TO_FILE, DataFrame)
    temp = sim_df[:,:CF]
    sim_data[:,i] = temp
end

# load real data

 _PATH_TO_REAL_DATA = joinpath(pwd(),"data")
 _PATH_TO_SORTED_REAL_DATA = joinpath(_PATH_TO_REAL_DATA, "sorted")

real_plot = CSV.read((joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-$(visit)-TEG-$(tpa)-nM-tPA.csv")),DataFrame)

# pull time vector out
real_time = Array{Float64}(real_plot[:,1])

# make it plottable
real_cf = Array{Float64}(real_plot[:,2:end])

# set path to figs
_PATH_TO_FIGS_VALIDATION = joinpath(pwd(),"figs_validation")
_PATH_TO_ENS_FIGS = joinpath(_PATH_TO_FIGS_VALIDATION,"ens")
path_to_CF_ens_figs_png = joinpath(_PATH_TO_ENS_FIGS, "tPA_$(tpa)nM_visit_$(visit)_CF_ens_validation.png") 
path_to_CF_ens_figs_pdf = joinpath(_PATH_TO_ENS_FIGS, "tPA_$(tpa)nM_visit_$(visit)_CF_ens_validation.pdf") 

# plot set of simulations for a given visit number, [tPA] combo
fig1 = plot(t, sim_data, xticks=0.0:15:90, yticks=0.0:15:90, xlim = (0,90), ylim = (0,90),label="",lw = 1.25,c=colorant"#89CCE2", bg="aliceblue", background_color_outside="white", framestyle = :box, xlabel="Time (min)", ylabel="CF (mm)")
       plot!(t,sim_data[:,1,], xticks=0.0:15:90, yticks=0.0:15:90, xlim = (0,90), ylim = (0,90),label="Simulation",lw = 1.25,c=colorant"#89CCE2", bg="aliceblue", background_color_outside="white", framestyle = :box, xlabel="Time (min)", ylabel="CF (mm)",fg_legend = :transparent)

       scatter!(real_time,real_cf, xticks=0.0:15:90,  xlim = (0,90), mc= :salmon, msc= :salmon, ma=0.3, ms=1,label="")
       scatter!(real_time,real_cf[:,1], xticks=0.0:15:90,  xlim = (0,90), mc= :salmon, msc= :salmon, ma=0.3, ms=1,label="Experimental")

savefig(fig1,path_to_CF_ens_figs_png)
savefig(fig1,path_to_CF_ens_figs_pdf)