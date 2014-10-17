require("phylo_model/phylo_mcmc.jl")
require("samplers/hmc.jl")
require("samplers/refractive_sampler.jl")

if !isdefined(:filename)
    println("Using default input file CLL077.csv")    
end
require("data_utils/read_phylosub_data.jl")

if !isdefined(:plotting)
    plotting=true
end

model_spec = ModelSpecification(ones(3)/3, false, false, plotting)

#model_spec.debug = true

lambda = 0.2
gamma = 1.0
alpha = 1.0

if filename == "CLL077.csv" || filename == "CLL003.csv"
    init_K = 4
elseif filename == "CLL006.csv"
    init_K = 5
end

if !isdefined(:num_trials)
    num_trials = 1
end

data = DataState(AA, DD, mu_r, mu_v, names)

result = mcmc(data, lambda, gamma, alpha, init_K, model_spec, 10000, 200) 

#results = run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, w_sigma, b_sigma, 1000, 500, num_trials, "prison_L$(lambda)_G$(gamma)_P$(trnpct)_PW$(positive_W)", "../results/prison/")

