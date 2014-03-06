require("model.jl")
require("data_utils.jl")
function train_test_split(YY::Array{Array{Float64, 2},1},
                          train_pct::Float64,
                          symmetric_split::Bool,
                          rand_seed::Int64)
    assert(train_pct <= 1.0 && train_pct >= 0.0)
    Y = YY[1] 
    Ytri = triu(Y)

    if symmetric_split
        y_inds = find(triu(ones(size(Y)),1))
    else
        y_inds = find(ones(size(Y))-eye(size(Y)[1]))
    end
    srand(rand_seed)
    shuffle!(y_inds)
    train_end = ifloor(train_pct * length(y_inds))

    y_train_inds = y_inds[1:train_end]

    train_mask = zeros(size(Y))
    train_mask[y_train_inds] = 1

    if symmetric_split
        train_mask = triu(train_mask)' + triu(train_mask,1)
    end
    test_mask = 1 - train_mask

    Ytrain = deepcopy(YY)
    Ytest = deepcopy(YY)

    for i = 1:length(YY)

        Ytrain[i][find(test_mask)] = -1
        Ytest[i][find(train_mask)] = -1

        # remove diagonal terms
        diag_mask = diagm(ones(size(Y)[1]))
        Ytrain[i][find(diag_mask)] = -1
        Ytest[i][find(diag_mask)] = -1
    end

    (Ytrain, Ytest)
end

function run_and_save(result_path, id_string, trial, mcmc_args...)
    results = mcmc(mcmc_args...)
    trial_string = "$(id_string)_$(trial)"
    models = results[end]
    model = models[end]
    data = mcmc_args[1]
    datas = Array(Any,2)

    # code currently assumes train/test splits are the same across all observations, so only save first observation
    datas[1] = data.Ytrain[1]
    datas[2] = data.Ytest[1]
    trees_array = [model2array(models[i])[1] for i = 1:length(models)]
    features_array = [model2array(models[i])[2] for i = 1:length(models)]
    weights_array = [model2array(models[i])[3] for i = 1:length(models)]
    save("$result_path/metrics_$trial_string.jla", results[1:end-1])
    save("$result_path/trees_$trial_string.jla", trees_array)
    save("$result_path/features_$trial_string.jla", features_array)
    save("$result_path/weights_$trial_string.jla", weights_array)
    save("$result_path/data_$trial_string.jla", datas)
end

function run_batch(model_spec::ModelSpecification,
                   Y::Array{Array{Float64,2},1},
                   symmetric_split::Bool,
                   train_pct::Float64,
                   lambda::Float64,
                   gamma::Float64,
                   w_sigma::Float64,
                   num_iterations::Int,
                   burnin_iterations::Int,
                   num_trials::Int,
                   id_string,
                   result_path)

    assert(burnin_iterations < num_iterations)

    N = size(Y[1],1)
    X_r = zeros((N,N,0))
    X_p = zeros((N,0))
    X_c = zeros((N,0))

    datas = Array(Any, num_trials) 
    for i = 1:num_trials
        Ytrain, Ytest = train_test_split(Y, train_pct, symmetric_split, i)
        datas[i] = DataState(Ytrain, Ytest, copy(X_r), copy(X_p), copy(X_c))
    end

    if num_trials > 1
        result_paths = [copy(result_path) for i = 1:num_trials]
        id_strings = [copy(id_string) for i = 1:num_trials]
        trials = [1:num_trials]

        lambdas = lambda*ones(num_trials)
        gammas = gamma*ones(num_trials)
        sigmas = w_sigma*ones(num_trials)
        model_specs = [deepcopy(model_spec) for i = 1:num_trials]
        num_iterses = num_iterations*ones(Int, num_trials)
        burn_iterses = burnin_iterations*ones(Int, num_trials)

        pmap(run_and_save, result_paths, id_strings, trials, datas, 
             lambdas, gammas, sigmas, model_specs, num_iterses, burn_iterses)

    else
        run_and_save(result_path, id_string, 1, datas[1], 
             lambda, gamma, w_sigma, model_spec, num_iterations, burnin_iterations)
    end


end

function restore_all(result_path, id_string, num_trials)
    trials = 1:num_trials 
    metrics = restore_in_range("$result_path/metrics_$id_string", trials)
    trees = restore_in_range("$result_path/trees_$id_string", trials)
    features = restore_in_range("$result_path/features_$id_string", trials)
    weights = restore_in_range("$result_path/weights_$id_string", trials)
    datas = restore_in_range("$result_path/data_$id_string", trials)

    results = Dict{ASCIIString, Any}()
    results["metrics"] = metrics_array2dict(metrics)
    results["trees"] = trees
    results["features"] = features
    results["weights"] = weights
    results["datas"] = datas

    results
end

function metrics_array2dict(metrics)
    dictionaries = {Dict{ASCIIString, Array}() for i = 1:length(metrics)}

    for i = 1:length(metrics)
        metric = metrics[i]
        dict = dictionaries[i]
        dict["iters"] = metric[1]
        dict["train_errors"] = metric[2]
        dict["test_errors"] = metric[3]
        dict["avg_test_LLs"] = metric[4]
        dict["AUCs"] = metric[5]
        dict["Ks"] = metric[6]
        dict["trainLLs"] = metric[7]
        dict["testLLs"] = metric[8]
    end
    dictionaries
end
