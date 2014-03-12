require("experiment_utils.jl")
require("mcmc.jl")
require("hmc.jl")
require("refractive_sampler.jl")

include("read_synthetic_data.jl")

opts = Dict{ASCIIString, Any}()
hmc_opts = @options numsteps=2 stepsize=0.01
opts["hmc"] = hmc_opts 
ref_opts = @options w=0.01 m=2 refractive_index_ratio=1.3

opts["RTJ_sampler"] = refractive_sampler
opts["RTJ_options"] = ref_opts
opts["W_sampler"] = hmc_sampler
opts["W_options"] = hmc_opts
#opts["RTJ_sampler"] = hmc_sampler
#opts["RTJ_options"] = hmc_opts

plotting = true
if !isdefined(:num_trials)
    num_trials=1
end

w_sigma = 100.0

model_spec = ModelSpecification(false, false, false, false, false, false, normal_logpdf, normal_logpdf_dx, 0.0, ones(3)/3, 1.0, 1.0, opts, false, false, plotting)
X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

b_sigma = 1.0

trnpct = 1.0
symmetric_split = false
lambda = 0.1
gamma = 1.5
run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, w_sigma, b_sigma, 500, 200, num_trials, "synthetic_L$(lambda)_G$(gamma)_P$(trnpct)", "../results/synthetic/")

