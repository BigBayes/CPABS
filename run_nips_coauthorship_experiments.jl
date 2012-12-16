load("mcmc.jl")
load("read_nips_data.jl")
load("profile.jl")

model_spec = ModelSpecification(false, false, false, true, 0.1, false)

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

mcmc(Y234, X_r, X_p, X_c, 0.000000002, 0.5, model_spec, 8)
