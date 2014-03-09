require("experiment_utils.jl")
require("mcmc.jl")
require("hmc.jl")
require("refractive_sampler.jl")

include("read_nips_data.jl")

tree_global_move_prob = 0.5
Z_sample_prob = 0.25

opts = Dict{ASCIIString, Any}()

positive_W = true
if !isdefined(:plotting)
    plotting = false
end

if positive_W
    hmc_opts = @options numsteps=4 stepsize=0.0002
    opts["hmc"] = hmc_opts
    ref_opts = @options w=0.2 m=2 refractive_index_ratio=1.3

    opts["RTJ_sampler"] = refractive_sampler
    opts["RTJ_options"] = ref_opts 
    model_spec = ModelSpecification(false, false, false, false, false, false, exp_logpdf, exp_logpdf_dx, log(0.01), ones(3)/3, tree_global_move_prob, Z_sample_prob, opts, false, false,plotting)
    model_spec.positive_W = true
else
    hmc_opts = @options numsteps=2 stepsize=0.0003
    opts["hmc"] = hmc_opts
    ref_opts = @options w=0.1 m=2 refractive_index_ratio=1.3

    opts["RTJ_sampler"] = refractive_sampler
    opts["RTJ_options"] = ref_opts 
    model_spec = ModelSpecification(false, false, false, false, false, false, normal_logpdf, normal_logpdf_dx, 0.0, ones(3)/3, tree_global_move_prob, Z_sample_prob, opts, false, false,plotting)
end
model_spec.diagonal_W = false

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

w_sigma = 1.0
b_sigma = 1.0

model_spec.symmetric_W = true
symmetric_split = true
model_spec.use_parenthood = true
model_spec.use_childhood = true

trnpct = 0.8
lambda = 1.0
gamma = 2.0
num_iterations = 500
run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, w_sigma, b_sigma, num_iterations, 200, 10, 
          "NIPS_D$(model_spec.diagonal_W)_S$(model_spec.positive_W)_L$(lambda)_G$(gamma)_T$(tree_global_move_prob)_P$(trnpct)_N$num_iterations", "../results/nips/")

# split into train/test
#Ytrain, Ytest = train_test_split(YY, .8)
#
#
#data = DataState(Ytrain, Ytest, X_r, X_p, X_c)
#
#results = mcmc(data, 0.01, 0.5, model_spec, 5)
#models = results[end]
#model = models[end]
## can't save models directly
#save("testfile.jlz", results[1:end-1])
#save("testmodel.jlz", model2array(model)
