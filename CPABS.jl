module CPABS

#using NumericExtensions
using Distributions
using Options
using OptionsMod
using HDF5
using JLD
using JSON

#import NumericExtensions: evaluate, result_type

include("model/phylo_mcmc.jl")
include("data_utils/read_phylosub_data.jl")
include("utils/parallel_utils.jl")
include("utils/general_macros.jl")

include("eval/eval_phylo.jl")
include("eval/phylo_plot.jl")
include("run_scripts/run_phylo_experiments.jl")

end
