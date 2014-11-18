require("phylo_model/phylo_mcmc.jl")
require("samplers/hmc.jl")
require("samplers/refractive_sampler.jl")
require("data_utils/read_phylosub_data.jl")

#if !isdefined(:filename)
#    println("Using default input file CLL077.csv")    
#    filename = "CLL077.csv"
#end

function run_all_emptysims_experiments(alpha)
    filenames = readdir("../data/phylosub/emptysims")

    for fname in filenames
        @spawn run_phylo_experiment("emptysims/$fname", alpha)
    end
end

function run_all_aldous_experiments(alpha)
    filenames = readdir("../data/phylosub/aldous")

    for fname in filenames
        @spawn run_phylo_experiment("aldous/$fname", alpha)
    end
end

function run_phylo_experiment(filename, alpha::Float64; 
                              wl_boundaries::Vector{Float64} = [Inf], #[-1400:50:-1250.0],
                              wl_K_boundaries::Vector{Float64} = [4,5,6,7,Inf],
                              wl_f0::Float64 = 1.0,
                              wl_histogram_test_ratio::Float64 = 0.3)

    (AA, DD, mu_r, mu_v, names) = read_phylosub_data(filename)

    (M, S) = size(AA)

    WL_state = WangLandauState(wl_boundaries, wl_K_boundaries, wl_f0, wl_histogram_test_ratio)

    if !isdefined(:plotting)
        model_spec = ModelSpecification(ones(3)/3, false, false, false, WL_state)
    else
        model_spec = ModelSpecification(ones(3)/3, false, false, plotting, WL_state)
    end


    #model_spec.debug = true

    lambda = 0.2
    gamma = 1.0

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
        jump_lag = 20
        jump_scan_length = 20
        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", filename)
        init_K = int(m.captures[1])-1
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
        rand_restarts=0
    elseif contains(filename, "aldous")
        filename_base="aldous"
        eta_Temp = 0.1
        jump_lag = 1
        jump_scan_length = 20
        m = match(r"\.([0-9]+)\.([0-9]+)\.", filename)
        init_K = 4
        D = int(m.captures[1])
        count = int(m.captures[2])
        rand_restarts=0
    end

    if !isdefined(:num_trials)
        num_trials = 1
    end

    data = DataState(AA, DD, mu_r, mu_v, names)

    result = mcmc(data, lambda, gamma, alpha, init_K, model_spec, 100000, 50, jump_lag = jump_lag, jump_scan_length = jump_scan_length, eta_Temp=eta_Temp, rand_restarts=rand_restarts)

    wl_state = model_spec.WL_state
    total_partition_mass = logsumexp(wl_state.partition_function)

    cocluster_matrix = zeros(M,M)
    (iters, Ks, trainLLs, models) = result
    sum_w = 0.0
    for n = 1:length(models)
        model = models[n]
        N::Int = (length(model.tree.nodes)+1)/2

        log_pdf = full_pdf(model, model_spec, data)
        w = get_partition_function(wl_state, N-1, log_pdf)
        w = exp(w - total_partition_mass)
        sum_w += w

        Z = model.Z
        C = [Int64[] for x=1:maximum(Z)]
        for i = 1:length(Z)
            push!(C[Z[i]], i)
        end
        for c in C
            cocluster_matrix[c,c] += w
        end
    end
    cocluster_matrix /= sum_w

    if filename_base == "emptysims"
        writedlm("../results/phylo/$filename_base.ccm.$alpha.$init_K.$D.$M_per_cluster.csv", cocluster_matrix, ",")
    elseif filename_base == "aldous"
        writedlm("../results/phylo/$filename_base.ccm.$alpha.$init_K.$D.$count.csv", cocluster_matrix, ",")
    else
        return cocluster_matrix
    end
end

