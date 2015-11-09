using HDF5
using JLD

@load "../data/synthetic/synthetic_data1.jld"

N = size(Y,1)

YY = Array(Array{Float64,2},1)
YY[1] = convert(Array{Float64},Y)
