function amplitude(time::Array{Float64,1},fibrin::Array{Float64,1},tPA::Float64,thrombin::Array{Float64,1},prothrombin::Float64)

    Rt = length(time)
    a0 = 0.125                                # base amplitude
    kx = 2500 - 1100*tPA                       # tPA function
    td = 4                                     # delay parameter in minutes

    CF = Array{Float64}(undef,(Rt,))            # clotting firmness
    

    for i ∈ 1:Rt
        fx = fibrin[i]
        FIIa = thrombin[i]

        if(time[i]>td)
            tau = 0.1*(1-(FIIa/prothrombin))
            S = 1.1+7.0*tPA                   # calculated using clot parameters
            a1 = S*(1 - exp(-tau*(time[i]-td)))
        else
            a1 = 0 
        end
    
        x = a0 + a1*(fx^2/(fx^2+kx^2))
        CF[i] = x

    end

    return CF
end