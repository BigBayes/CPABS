require("experiment_utils.jl")
require("mcmc.jl")
include("read_monks_data.jl")

model_spec = ModelSpecification(false, false, false, false, false, ones(3)/3, 1.0, 1.0, false, false)
X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

trnpct = 0.8
symmetric_split = false
lambda = 0.01
gamma = 0.9
run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, 1.0, 500, 200, 10, "monks_L$(lambda)_G$(gamma)_P$(trnpct)", "../results/monks/")

