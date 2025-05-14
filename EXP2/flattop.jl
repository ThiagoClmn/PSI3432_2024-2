function flattop(N::Int)
    n = 0:N-1
    L = (N-1)/2
    w = 0.24726 .+ 0.46071cos.(2pi*(n.-L)/(N-1)) .+ 0.25078cos.(4pi*(n.-L)/(N-1)) .+ 0.04125cos.(6pi*(n.-L)/(N-1))
    return w 
end