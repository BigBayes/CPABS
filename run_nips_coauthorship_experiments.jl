require("experiment_utils.jl")
require("mcmc.jl")
include("read_nips_data.jl")

tree_global_move_prob = 0.5
model_spec = ModelSpecification(false, false, false, false, false, ones(3)/3, tree_global_move_prob, 1.0, false, false)


X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

symmetric_split = true
trnpct = 0.8
lambda = 0.01
gamma = 0.95
num_iterations = 500
run_batch(model_spec, YY, symmetric_split, trnpct, lambda, gamma, 1.0, num_iterations, 200, 10, 
          "NIPS_L$(lambda)_G$(gamma)_T$(tree_global_move_prob)_P$(trnpct)_N$num_iterations", "../results/nips/")

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
