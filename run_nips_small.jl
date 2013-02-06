@everywhere require("mcmc.jl")
@everywhere require("data_utils.jl")
@everywhere require("read_nips_data.jl")
@everywhere require("experiment_utils.jl")

model_spec = ModelSpecification(false, false, false, false, false, ones(3)/3, 0.2, 1.0, false, false)

model_spec.symmetric_W = false
model_spec.diagonal_W = false
model_spec.debug = false

X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

perm = randperm(234)
permute_rows_and_cols!(YY,perm)

Ysmall = YY[1:70,1:70]

trnpct = 0.8
symmetric_split = true
run_batch(model_spec, Ysmall, symmetric_split, trnpct, 0.3, 0.5, 300, 150, 1, "NIPS_small_$trnpct", "../results/nips_small/")

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
