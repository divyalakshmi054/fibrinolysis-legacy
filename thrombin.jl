function thrombin(time::Float64, FIIa_itp)
    #include -
    return FIIa_itp(time*100+1)
end