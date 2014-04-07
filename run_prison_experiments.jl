require("experiment_utils.jl")
require("mcmc.jl")
require("hmc.jl")
require("refractive_sampler.jl")

include("read_prison_data.jl")

opts = Dict{ASCIIString, Any}()
hmc_opts = @options numsteps=4 stepsize=0.05
opts["hmc"] = hmc_opts 
ref_opts = @options w=0.05 m=2 refractive_index_ratio=1.3

opts["RTJ_sampler"] = refractive_sampler
opts["RTJ_options"] = ref_opts
opts["W_sampler"] = hmc_sampler
opts["W_options"] = hmc_opts
opts["ab_sampler"] = hmc_sampler
opts["ab_options"] = hmc_opts

positive_W = false
if !isdefined(:plotting)
    plotting = false
end

model_spec = ModelSpecification(false, false, false, false, false, false, normal_logpdf, normal_logpdf_dx, 0.0, ones(3)/3, 1.0, 1.0, opts, false, false, plotting)

model_spec.use_parenthood = true
model_spec.use_childhood = true
#model_spec.debug = true

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

w_sigma = 1.0
b_sigma = 0.1

if !isdefined(:trnpct)
    trnpct = 0.8
end
symmetric_split = false
lambda = 0.1
gamma = 2.0

if !isdefined(:num_trials)
    num_trials = 10
end

run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, w_sigma, b_sigma, 1000, 500, num_trials, "prison_L$(lambda)_G$(gamma)_P$(trnpct)_PW$(positive_W)", "../results/prison/")

