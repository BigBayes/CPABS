require("phylo_model/phylo_mcmc.jl")
require("samplers/hmc.jl")
require("samplers/refractive_sampler.jl")
require("data_utils/read_phylosub_data.jl")

#if !isdefined(:filename)
#    println("Using default input file CLL077.csv")    
#    filename = "CLL077.csv"
#end

function run_all_emptysims_experiments()
    filenames = readdir("../data/phylosub/emptysims")

    for fname in filenames
        @spawn run_phylo_experiment(fname)
    end
end


function run_phylo_experiment(filename)
    (AA, DD, mu_r, mu_v, names) = read_phylosub_data(filename)

    (M, S) = size(AA)
    if !isdefined(:plotting)
        plotting=false
    end

    model_spec = ModelSpecification(ones(3)/3, false, false, plotting)

    #model_spec.debug = true

    lambda = 0.2
    gamma = 1.0
    alpha = 1.0

    #if filename == "CLL077.csv" || filename == "CLL003.csv"
    #    init_K = 4
    #elseif filename == "CLL006.csv"
    #    init_K = 5
    #end

    if contains(filename, "CLL")
        init_K=4
        filename_base="CLL"
        eta_Temp = 0.0000001
        jump_lag = Inf #100
        jump_scan_length = 100
    elseif contains(filename, "emptysims")
        filename_base="emptysims"
        eta_Temp = 0.001
        jump_lag = Inf
        jump_scan_length = 20
        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", filename)
        init_K = int(m.captures[1])-1
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
    end

    if !isdefined(:num_trials)
        num_trials = 1
    end

    data = DataState(AA, DD, mu_r, mu_v, names)

    result = mcmc(data, lambda, gamma, alpha, init_K, model_spec, 10000, 500, jump_lag = jump_lag, jump_scan_length = jump_scan_length, eta_Temp=eta_Temp, rand_restarts=10)

    cocluster_matrix = zeros(M,M)
    (iters, Ks, trainLLs, models) = result
    for n = 1:length(models)
        model = models[n]
        Z = model.Z
        C = [Int64[] for x=1:maximum(Z)]
        for i = 1:length(Z)
            push!(C[Z[i]], i)
        end
        for c in C
            cocluster_matrix[c,c] += 1
        end
    end
    cocluster_matrix /= length(models)

    if filename_base == "emptysims" || filename_base == "aldous"
        writedlm("../results/phylo/$filename_base.ccm.$init_K.$D.$M_per_cluster.csv", cocluster_matrix, ",")
    else
        return cocluster_matrix
    end
end

