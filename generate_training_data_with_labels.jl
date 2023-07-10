# Where are we?
_BASE_PATH = pwd()
_PATH_TO_DATA = joinpath(_BASE_PATH, "data")

# Package that we need -
using DataFrames
using CSV
using Distributions
using Statistics
using LinearAlgebra

# load the actual data -
path_to_experimental_data = joinpath(_PATH_TO_DATA, "Training-Thrombin-TF-Hormones-Fibrinolysis.csv")
full_df = CSV.read(path_to_experimental_data, DataFrame)
(number_of_records, _) = size(full_df) 

# what are the col names? (exclude subject id and visit id)
data_col_symbols = Symbol.(names(full_df)[1:end])
number_of_fields = length(data_col_symbols)

# setup scale factor dictionary to convert to concentration units -
SF = 1e9
scale_factor_dict = Dict()
scale_factor_dict[:II] = (1.4e-6)*SF*(1/100)
scale_factor_dict[:V] = (2e-8)*SF*(1/100)
scale_factor_dict[:VII] = (1e-8)*SF*(1/100)
scale_factor_dict[:VIII] = (7e-10)*SF*(1/100)
scale_factor_dict[:IX] = (9e-8)*SF*(1/100)
scale_factor_dict[:X] = (1.6e-7)*SF*(1/100)
scale_factor_dict[:XI] = (1e-8)*SF*(1/100)
scale_factor_dict[:XII] = (1e-8)*SF*(1/100)
scale_factor_dict[:AT] = (3.4e-6)*SF*(1/100)
scale_factor_dict[:PC] = (6.3e-8)*SF*(1/100)
scale_factor_dict[:TFPI] = 1.0
scale_factor_dict[:Lagtime] = 1.0
scale_factor_dict[:Peak] = 1.0
scale_factor_dict[:TPeak] = 1.0
scale_factor_dict[:Max] = 1.0
scale_factor_dict[:EPT] = 1.0
scale_factor_dict[:Fbgn] = 1.0  #calculated using molecular mass of compound
scale_factor_dict[:Plgn] = (2e-6)*SF*(1/100) #calculated using molecular mass of compound
scale_factor_dict[:PAI1] = 1.0
scale_factor_dict[:TAFI] = 75.0*(1/100)
scale_factor_dict[:XIII] = 68.0*(1/100)
# scale_factor_dict[:a2AP] = (1.53e-8)*SF #calculated using molecular mass of compound - let's revisit this later
scale_factor_dict[:subjid] = 1
scale_factor_dict[:visitid] = 1
scale_factor_dict[:CT] = 1.0
scale_factor_dict[:MCF] = 1.0
scale_factor_dict[:alpha] = 1.0
scale_factor_dict[:LOT] = 1.0
scale_factor_dict[:ML] = 1.0
scale_factor_dict[:LT] = 1.0
scale_factor_dict[:AUC] = 1.0

# initialize -
transformed_df = DataFrame()
for rᵢ ∈ 1:number_of_records


    tmp = Vector()
    for fᵢ ∈ 1:number_of_fields
        field_symbol = data_col_symbols[fᵢ]
        value = scale_factor_dict[field_symbol]*full_df[rᵢ,field_symbol]
        push!(tmp, value)
    end

    # create new record -
    transformed_tuple = (
        subjid = tmp[1], 
        visitid = tmp[2],
        II = tmp[3],
        V = tmp[4],
        VII = tmp[5],
        VIII = tmp[6],
        IX = tmp[7],
        X = tmp[8],
        XI = tmp[9],
        XII = tmp[10],
        AT = tmp[11],
        PC = tmp[12],
        TFPI = tmp[13],
        Lagtime = tmp[14],
        Peak = tmp[15],
        TPeak = tmp[16],
        Max = tmp[17],
        EPT = tmp[18],
        Fbgn = tmp[19],
        Plgn = tmp[20],
        PAI1 = tmp[21],
        TAFI = tmp[22],
        FXIII = tmp[23],
        CT = tmp[24],
        MCF = tmp[25],
        alpha = tmp[26],
        LOT = tmp[27],
        ML = tmp[28],
        LT = tmp[29],
        AUC = tmp[30]
    )
    push!(transformed_df, transformed_tuple)
end

# dump sample data to disk -
path_to_transformed_data = joinpath(_PATH_TO_DATA, "Training-Fibrinolysis-Transformed-w-Labels.csv")
CSV.write(path_to_transformed_data, transformed_df)