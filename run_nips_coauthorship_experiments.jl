load("mcmc.jl")
load("read_nips_data.jl")
load("profile.jl")

model_spec = ModelSpecification(false, false, false, true, false, ones(3)/3, 0.2, 1.0, true, true)
X_r = zeros((0,0,0))
X_p = zeros((0,0))
X_c = zeros((0,0))

trainmask = zeros(Int64, size(Y234))

# split into train/test
y_inds = linspace(1,length(Y234),length(Y234))
shuffle!(y_inds)
train_end = ifloor(.8 * length(Y234))

y_train_inds = y_inds[1:train_end]
y_test_inds = y_inds[train_end+1:end]

y_ind_nnz = find(Y234) 
Ytrain = copy(Y234)
Ytrain[y_ind_nnz] = 1
Ytest = copy(Ytrain)

Ytrain[y_test_inds] = -1
Ytest[y_train_inds] = -1


data = DataState(Ytrain, Ytest, X_r, X_p, X_c)

models = mcmc(data, 0.01, 0.5, model_spec, 500)
