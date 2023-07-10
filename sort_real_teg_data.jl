include("Include.jl")

 _PATH_TO_REAL_DATA = joinpath(pwd(),"data")
real_df = CSV.read(joinpath(_PATH_TO_REAL_DATA,"Training-TEG-Transpose.csv"),DataFrame)

# clean up missing => for now
real_df = mapcols(col -> replace(col, missing => 0), real_df)

# filter by yes or no tPA

no_tpa_df = filter(:TPA => x->(x=="N"), real_df)
yes_tpa_df = filter(:TPA => x->(x=="Y"), real_df)

# filter by visit

v1_no_df = filter(:Visit => x->(x=="V1"), no_tpa_df) 
v2_no_df = filter(:Visit => x->(x=="V2"), no_tpa_df) 
v4_no_df = filter(:Visit => x->(x=="V4"), no_tpa_df) 
v1_yes_df = filter(:Visit => x->(x=="V1"), yes_tpa_df) 
v2_yes_df = filter(:Visit => x->(x=="V2"), yes_tpa_df) 
v4_yes_df = filter(:Visit => x->(x=="V4"), yes_tpa_df) 

# get real data in plottable form

# convert given time(s) to min
real_t = collect(5:5:10800)/60

# rearranging all visits, no tPA

no_tpa_df.Time_sec = string.(1:nrow(no_tpa_df))
no_tpa_plot = permutedims(no_tpa_df[:,4:end],"Time_sec")
no_tpa_plot = no_tpa_plot[:,2:end]
no_tpa_plot = insertcols!(no_tpa_plot, 1, :Time_sec => real_t)

# rearranging all visits, with tPA

yes_tpa_df.dot = string.(1:nrow(yes_tpa_df))
yes_tpa_plot = permutedims(yes_tpa_df[:,4:end],"dot")
yes_tpa_plot = yes_tpa_plot[:,2:end]
yes_tpa_plot = insertcols!(yes_tpa_plot, 1, :Time_sec => real_t)


# rearranging visit 1 dataframes

v1_no_df.Time_sec = string.(1:nrow(v1_no_df))
v1_no_tpa_plot = permutedims(v1_no_df[:,4:end],"Time_sec")
v1_no_tpa_plot = v1_no_tpa_plot[:,2:end]
v1_no_tpa_plot = insertcols!(v1_no_tpa_plot, 1, :Time_sec => real_t)


v1_yes_df.Time_sec = string.(1:nrow(v1_yes_df))
v1_yes_tpa_plot = permutedims(v1_yes_df[:,4:end],"Time_sec")
v1_yes_tpa_plot = v1_yes_tpa_plot[:,2:end]
v1_yes_tpa_plot = insertcols!(v1_yes_tpa_plot, 1, :Time_sec => real_t)

# rearranging visit 2 dataframes

v2_no_df.Time_sec = string.(1:nrow(v2_no_df))
v2_no_tpa_plot = permutedims(v2_no_df[:,4:end],"Time_sec")
v2_no_tpa_plot = v2_no_tpa_plot[:,2:end]
v2_no_tpa_plot = insertcols!(v2_no_tpa_plot, 1, :Time_sec => real_t)

v2_yes_df.Time_sec = string.(1:nrow(v2_yes_df))
v2_yes_tpa_plot = permutedims(v2_yes_df[:,4:end],"Time_sec")
v2_yes_tpa_plot = v2_yes_tpa_plot[:,2:end]
v2_yes_tpa_plot = insertcols!(v2_yes_tpa_plot, 1, :Time_sec => real_t)

# rearranging visit 4 dataframes

v4_no_df.Time_sec = string.(1:nrow(v4_no_df))
v4_no_tpa_plot = permutedims(v4_no_df[:,4:end],"Time_sec")
v4_no_tpa_plot = v4_no_tpa_plot[:,2:end]
v4_no_tpa_plot = insertcols!(v4_no_tpa_plot, 1, :Time_sec => real_t)

v4_yes_df.Time_sec = string.(1:nrow(v4_yes_df))
v4_yes_tpa_plot = permutedims(v4_yes_df[:,4:end],"Time_sec")
v4_yes_tpa_plot = v4_yes_tpa_plot[:,2:end]
v4_yes_tpa_plot = insertcols!(v4_yes_tpa_plot, 1, :Time_sec => real_t)

# dump sorted real data
_PATH_TO_SORTED_REAL_DATA = joinpath(_PATH_TO_REAL_DATA, "sorted")

CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-allvisit-TEG-0-nM-tPA.csv"), no_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-allvisit-TEG-4-nM-tPA.csv"), yes_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-1-TEG-0-nM-tPA.csv"), v1_no_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-1-TEG-4-nM-tPA.csv"), v1_yes_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-2-TEG-0-nM-tPA.csv"), v2_no_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-2-TEG-4-nM-tPA.csv"), v2_yes_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-4-TEG-0-nM-tPA.csv"), v4_no_tpa_plot)
CSV.write(joinpath(_PATH_TO_SORTED_REAL_DATA, "REAL-visit-4-TEG-4-nM-tPA.csv"), v4_yes_tpa_plot)