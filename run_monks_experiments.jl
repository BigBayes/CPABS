require("experiment_utils.jl")
require("mcmc.jl")
include("read_monks_data.jl")

model_spec = ModelSpecification(false, false, false, false, false, ones(3)/3, 1.0, 1.0, false, false)
X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

trnpct = 0.8
symmetric_split = false
run_batch(model_spec, YY, symmetric_split, trnpct, 0.5, 0.5, 1.0, 100, 50, 1, "monks_$trnpct", "../results/monks/")
