@everywhere load("mcmc.jl")
@everywhere load("data_utils.jl")
load("read_nips_data.jl")
load("profile.jl")
@everywhere load("experiment_utils.jl")

model_spec = ModelSpecification(false, false, false, true, false, ones(3)/3, 0.2, 1.0, false, false)
X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

run_batch(model_spec, YY, 0.8, 0.5, 0.5, 500, 200, 10, "NIPS_test", "../results")

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
