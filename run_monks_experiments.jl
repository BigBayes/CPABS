require("experiment_utils.jl")
require("mcmc.jl")
require("hmc.jl")
require("refractive_sampler.jl")

include("read_monks_data.jl")

w_sigma = 1.0

opts = Dict{ASCIIString, Any}()
hmc_opts = @options numsteps=4 stepsize=0.02
opts["hmc"] = hmc_opts 
ref_opts = @options w=0.2 m=1 refractive_index_ratio=1.3

#opts["RTJ_sampler"] = refractive_sampler
#opts["RTJ_options"] = ref_opts
opts["RTJ_sampler"] = hmc_sampler
opts["RTJ_options"] = hmc_opts


model_spec = ModelSpecification(false, false, false, false, false, false, ()->rand(Normal(0.0,w_sigma)), normal_logpdf, normal_logpdf_dx, ones(3)/3, 1.0, 1.0, opts, false, false)

model_spec.use_parenthood = true
model_spec.use_childhood = true

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

trnpct = 0.8
symmetric_split = false
lambda = 0.1
gamma = 2.0

if !isdefined(:num_trials)
    num_trials = 10
end

run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, w_sigma, 500, 200, num_trials, "monks_L$(lambda)_G$(gamma)_P$(trnpct)", "../results/monks/")

