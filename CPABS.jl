module CPABS

using ArgParse
using NumericExtensions
using Distributions
require("Options")
using OptionsMod
using HDF5, JLD


import NumericExtensions: evaluate, result_type

include("run_scripts/run_phylo_experiments.jl")



end
