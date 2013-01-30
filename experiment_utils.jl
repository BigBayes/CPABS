load("mcmc.jl")
load("read_nips_data.jl")
load("data_utils.jl")

function train_test_split(Y::Array{Int64, 2},
                          train_pct::Float64)
    assert(train_pct <= 1.0 && train_pct >= 0.0)
    
    y_inds = linspace(1,length(Y),length(Y))
    shuffle!(y_inds)
    train_end = ifloor(train_pct * length(Y))

    y_train_inds = y_inds[1:train_end]
    y_test_inds = y_inds[train_end+1:end]

    Ytrain = copy(Y)
    Ytest = copy(Y)
    Ytrain[y_test_inds] = -1
    Ytest[y_train_inds] = -1

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
    Ytrain, Ytest = train_test_split(Y, train_pct)
    X_r = zeros((0,0,0))
    X_p = zeros((0,0))
    X_c = zeros((0,0))

    data = DataState(Ytrain, Ytest, X_r, X_p, X_c)

    result_refs = Array(Any, num_trials)
    results = Array(Any, num_trials)

    result_paths = [copy(result_path) for i = 1:num_trials]
    id_strings = [copy(id_string) for i = 1:num_trials]
    train_pcts = train_pct * ones(num_trials)
    trials = linspace(1,num_trials,num_trials)

    datas = [copy(data) for i = 1:num_trials]
    lambdas = lambda*ones(num_trials)
    gammas = gamma*ones(num_trials)
    model_specs = [copy(model_spec) for i = 1:num_trials]
    num_iterses = num_iterations*ones(Int, num_trials)
    burn_iterses = burnin_iterations*ones(Int, num_trials)

    pmap(run_and_save, result_paths, id_strings, train_pcts, trials, datas, 
         lambdas, gammas, model_specs, num_iterses, burn_iterses)

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

function run_and_save(result_path, id_string, train_pct, trial, mcmc_args...)
    results = mcmc(mcmc_args...)
    trial_string = "$(id_string)_$(train_pct)_$(trial)"
    models = results[end]
    model = models[end]
    save("$result_path/metrics_$trial_string.jla", results[1:end-1])
    save("$result_path/model_$trial_string.jla", model2array(model))
end
