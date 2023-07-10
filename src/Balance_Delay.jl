# include -
include("Include.jl")

function _balances_delay_term(dx, x, p, t)

    # delay term -
    # delay_term = 1 - 1/(1+exp(0.1*(t-60)))

    # grab data from the parameter vector 
    α = p[1]                            # rate constant vector                        
    # α[1] = delay_term                 # FIIa
    G = p[2]                            # exponent array
    S = p[3]                            # stoichoimetric array
    number_of_dynamic_states = p[4]     # number of dynamic states
    static_factors_array = p[5]         # list of static factors 
    u = p[6]                            # get the input function

    # check: do we have a callback function?
    if (u === nothing)
        u = _callback;
    end
        
    # build the "state" array (dynamic | static)
    static_factors_array[5] = thrombin(t,FIIa_itp)
    state_array = vcat(x,static_factors_array)

    # compute the kinetics - powerlaw
    rV = powerlaw(state_array,α,G)

    # compute the rhs -> store in a temp vector
    tmp = S*rV + u(t,x,p)

    # populate the dx vector with the tmp vector -
    for i ∈ (1:number_of_dynamic_states)
        dx[i] = tmp[i]
    end    
end


function _callback(t::Float64, x::Array{Float64,1}, p::Array{Any,1})::Array{Float64,1}
    return zeros(length(x));
end