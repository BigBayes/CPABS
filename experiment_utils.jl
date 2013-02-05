require("mcmc.jl")
require("data_utils.jl")

function train_test_split(Y::Array{Int64, 2},
                          train_pct::Float64,
                          symmetric_split::Bool)
    assert(train_pct <= 1.0 && train_pct >= 0.0)
   
    Ytri = triu(Y)

    if symmetric_split
        y_inds = find(triu(ones(size(Y)),1))
    else
        y_inds = find(ones(size(Y))-eye(size(Y)[1]))
    end
    shuffle!(y_inds)
    train_end = ifloor(train_pct * length(y_inds))

    y_train_inds = y_inds[1:train_end]

    train_mask = zeros(size(Y))
    train_mask[y_train_inds] = 1

    if symmetric_split
        train_mask = triu(train_mask)' + triu(train_mask,1)
    end
    test_mask = 1 - train_mask

    Ytrain = copy(Y)
    Ytest = copy(Y)
    Ytrain[find(test_mask)] = -1
    Ytest[find(train_mask)] = -1

    # remove diagonal terms
    diag_mask = diagm(ones(size(Y)[1]))
    Ytrain[find(diag_mask)] = -1
    Ytest[find(diag_mask)] = -1

    (Ytrain, Ytest)
end

function run_batch(model_spec::ModelSpecification,
                   Y::Array{Int64,2},
                   train_pct::Float64,
                   lambda::Float64,
                   gamma::Float64,
                   num_iterations::Int,
                   burnin_iterations::Int,
                   num_trials::Int,
                   id_string,
                   result_path)

    assert(burnin_iterations < num_iterations)

    X_r = zeros((0,0,0))
    X_p = zeros((0,0))
    X_c = zeros((0,0))

    datas = Array(Any, num_trials) 
    for i = 1:num_trials
        Ytrain, Ytest = train_test_split(Y, train_pct, model_spec.symmetric_W)
        datas[i] = DataState(Ytrain, Ytest, copy(X_r), copy(X_p), copy(X_c))
    end

    result_refs = Array(Any, num_trials)
    results = Array(Any, num_trials)

    result_paths = [copy(result_path) for i = 1:num_trials]
    id_strings = [copy(id_string) for i = 1:num_trials]
    trials = linspace(1,num_trials,num_trials)

    lambdas = lambda*ones(num_trials)
    gammas = gamma*ones(num_trials)
    model_specs = [copy(model_spec) for i = 1:num_trials]
    num_iterses = num_iterations*ones(Int, num_trials)
    burn_iterses = burnin_iterations*ones(Int, num_trials)

    pmap(run_and_save, result_paths, id_strings, trials, datas, 
         lambdas, gammas, model_specs, num_iterses, burn_iterses)

#    run_and_save(result_paths[1], id_strings[1], train_pcts[1], trials[1], datas[1],
#         lambdas[1], gammas[1], model_specs[1], num_iterses[1], burn_iterses[1])

#    for trial = 1:num_trials
#        result_refs[trial] = @spawn mcmc(data, lambda, gamma, model_spec,
#                                         num_iterations, burnin_iterations)
#        println(result_refs[trial]) 
#    end

#    @parallel for trial = 1:num_trials
#        run_and_save(result_path, id_string, train_pct, trial,
#                   datas[trial], lambda, gamma, model_spec, num_iterations, burnin_iterations)
#    end
#    for trial = 1:num_trials
#        trial_string = "$(id_string)_$(train_pct)_$(trial)"
#        results[trial] = fetch(result_refs[trial])
#        println(results[trial][1:end-1])
#        models = results[trial][end]
#        model = models[end]
#        save("$result_path/metrics_$trial_string.jla", results[trial][1:end-1])
#        save("$result_path/model_$trial_string.jla", model2array(model))
#    end

end

function run_and_save(result_path, id_string, trial, mcmc_args...)
    results = mcmc(mcmc_args...)
    trial_string = "$(id_string)_$(trial)"
    models = results[end]
    model = models[end]
    data = mcmc_args[1]
    datas = Array(Any,2)
    datas[1] = data.Ytrain
    datas[2] = data.Ytest
    trees_array = [model2array(models[i])[1] for i = 1:length(models)]
    features_array = [model2array(models[i])[2] for i = 1:length(models)]
    weights_array = [model2array(models[i])[3] for i = 1:length(models)]
    save("$result_path/metrics_$trial_string.jla", results[1:end-1])
    save("$result_path/trees_$trial_string.jla", trees_array)
    save("$result_path/features_$trial_string.jla", features_array)
    save("$result_path/weights_$trial_string.jla", weights_array)
    save("$result_path/data_$trial_string.jla", datas)
end
