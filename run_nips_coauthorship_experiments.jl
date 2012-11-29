load("mcmc.jl")
load("read_nips_data.jl")
load("profile.jl")

model_spec = ModelSpecification(false, false, false, true)

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

mcmc(Y234, X_r, X_p, X_c, 0.2, 0.5, model_spec, 4)
