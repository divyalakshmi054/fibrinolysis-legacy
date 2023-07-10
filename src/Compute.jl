# let's get some properties of the CF v t curve to help with our training -

# load the include -
include("Include.jl")

function clot_properties()
    # load the training data -
    _PATH_TO_DATA = joinpath(pwd(),"data")

    # specify visit & tPA -
    visit = 4;  # visit number
    tPA = 4;    # tPA conc - should be an integer (0 or 4)
    path_to_training_data = joinpath(_PATH_TO_DATA, "sorted","REAL-visit-$(visit)-TEG-$(tPA)-nM-tPA.csv")
    unclean_training_df = CSV.read(path_to_training_data, DataFrame, copycols = true)

    # dimensions -
    (R,C) = size(unclean_training_df)
    training_df = Array{Float64}(undef,(R,C))

    # some data is 'missing' - let's replace that with 0 for now
    for i ∈ 1:R
        for j ∈ 1:C
            if(ismissing(unclean_training_df[i,j]))
                training_df[i,j] = 0
            else
                training_df[i,j] = unclean_training_df[i,j]
            end
        end
    end

    T = training_df[:,1]                            # define time
    data_vector = Array{Float64}(undef,(5,(C-1)))   # 5 properties, C-1 patients

    for i ∈ 1:(C-1)                                 # omitting column 1: time
        X = training_df[1:end,i+1]
        data_vector[:,i] = model_output_vector(T,X)
        print(data_vector[:,i],"\n")
    end
    print(size(data_vector))
    #time to dump, finally! -
    data_output_header = ["CT", "CFT", "MCF","MCFt", "AUC"]
    CSV.write(joinpath(_PATH_TO_DATA,"sorted","clot_parameters","Training-Clot-Parameters-REAL-visit-$(visit)-$(tPA)-nM-tPA.csv"),Tables.table(transpose(data_vector));header = data_output_header)
end

function model_output_vector(T::Array{Float64,1},X::Array{Float64,1})::Array{Float64,1}
    output_vector = Array{Float64,1}()

    #properties -
    CT = clot_time(x->x>=2.0, T, X)
    MCF = maximum(X)
    LOT = lot(MCF,T,X)
    LT = lysis_time(MCF,T,X)
    alpha_slope = compute_alpha_slope(MCF, T, X)
    area_under_curve = auc(T,X)

    #it's go time -
    push!(output_vector,CT)
    push!(output_vector,MCF)
   # push!(output_vector,LOT)
   # push!(output_vector,LT)
    push!(output_vector,alpha_slope)
    push!(output_vector,area_under_curve)
    return output_vector
end

function auc(T::Array{Float64,1}, X::Array{Float64,1})::Float64
    return integrate(T,X)
end

function lot(MCF::Float64,T::Array{Float64,1}, X::Array{Float64,1})::Float64
    idx_MCF = findfirst(x->x==MCF,X)
    idx_LOT = findfirst(x->x<=0.85*MCF,X[idx_MCF:end])
    if(isnothing(idx_LOT) == true)
        return 180.0
    end
    return (T[idx_LOT]+T[idx_MCF])
end

function lysis_time(MCF::Float64,T::Array{Float64,1},X::Array{Float64,1})::Float64
    idx_MCF = findfirst(x->x==MCF,X)
    idx_LT = findfirst(x->x<=0.10*MCF,X[idx_MCF:end])
    if(isnothing(idx_LT) == true)
        return 180.0
    end
    return (T[idx_LT]+T[idx_MCF])
end

function clot_time(rule::Function, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function clot_formation_time(rule::Function, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function compute_alpha_slope(MCF::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64
    idx_CT = findfirst(x->x>=2.0,X)
    if(MCF>=20.0)
        idx_CFT = findfirst(x->x>=20.0,X)
    elseif(MCF<20.0)
        idx_CFT = findfirst(x->x==MCF,X)
    end

    if(isnothing(idx_CT) == true || isnothing(idx_CFT) == true)
        return 60.0
    end
    slope = (X[idx_CFT] - X[idx_CT])/(T[idx_CFT] - T[idx_CT])
    alpha_angle = atan(slope)*(180/pi)
    return alpha_angle
end

function max_lysis(T::Array{Float64,1},X::Array{Float64,1})
    MCF = maximum(X)
    idx_90 = findfirst(x->x==90.0,T)
    ML = ((MCF-X[idx_90])/MCF)*100.0
    return ML
end

function compute_amp_30(CT::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(x->x>=(CT+30),T)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return X[idx]
end

function MCF_time(MCF::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(x->x==MCF,X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end