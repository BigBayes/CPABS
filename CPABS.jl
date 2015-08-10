module CPABS

using ArgParse
using NumericExtensions
using Distributions
require("Options")
using OptionsMod
using HDF5, JLD


import NumericExtensions: evaluate, result_type

include("tree.jl")
include("wang_landau.jl")

include("run_scripts/run_phylo_cmd.jl")

end
