# include -
include("Include.jl")

# PRIVATE BELOW

function _evaluate_w_delay(model::Dict{String,Any}; 
    tspan::Tuple{Float64,Float64} = (0.0,20.0), Δt::Float64 = 0.01, input::Union{Nothing,Function} = nothing)

    # get stuff from model -
    xₒ = model["initial_condition_array"]

    # build parameter vector -
    p = Array{Any,1}(undef,6)
    p[1] = model["α"]
    p[2] = model["G"]
    p[3] = model["S"]
    p[4] = model["number_of_dynamic_states"]
    p[5] = model["static_factors_array"]
    p[6] = input;
    # setup the solver -
    prob = ODEProblem(_balances_delay_term, xₒ, tspan, p; saveat = Δt)
    soln = solve(prob)

    # get the results from the solver -
    T = soln.t
    tmp = soln.u
    # build soln array -
    number_of_time_steps = length(T)
    number_of_dynamic_states = model["number_of_dynamic_states"]
    X = Array{Float64,2}(undef, number_of_time_steps,  number_of_dynamic_states);

    for i ∈ 1:number_of_time_steps
        soln_vector = tmp[i]
        for j ∈ 1:number_of_dynamic_states
            X[i,j] = soln_vector[j]
        end
    end

    # return -
    return (T,X)
end
# PRIVATE ABOVE

# PUBLIC BELOW
function evaluate_w_delay(model::BSTModel; tspan::Tuple{Float64,Float64} = (0.0,20.0), Δt::Float64 = 0.01, 
    input::Union{Nothing,Function} = nothing)::Tuple{Array{Float64,1}, Array{Float64,2}}

    try
        # convert the model object to the internal_model_dictionary -
        internal_model_dictionary = _build_BST(model);
        return _evaluate_w_delay(internal_model_dictionary, tspan = tspan, Δt = Δt, input = input); 

    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # throw -
        throw(vl_error_obj);
    end
end

function _build_BST(model::BSTModel)::Dict{String,Any}

    # initialize -
    internal_model_dictionary = Dict{String,Any}();

    # get the data from the model object, put into dictionary -
    internal_model_dictionary["number_of_dynamic_states"] = model.number_of_dynamic_states
    internal_model_dictionary["number_of_static_states"] = model.number_of_static_states
    internal_model_dictionary["list_of_dynamic_species"] = model.list_of_dynamic_species
    internal_model_dictionary["list_of_static_fators"] = model.list_of_static_species
    internal_model_dictionary["total_species_list"] = model.total_species_list
    internal_model_dictionary["static_factors_array"] = model.static_factors_array
    internal_model_dictionary["initial_condition_array"] = model.initial_condition_array
    internal_model_dictionary["list_of_reactions"] = model.list_of_reactions
    internal_model_dictionary["S"] = model.S
    internal_model_dictionary["G"] = model.G
    internal_model_dictionary["α"] = model.α

    # return -
    return internal_model_dictionary;
end