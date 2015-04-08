require("phylo_model/phylo_mcmc.jl")
require("samplers/hmc.jl")
require("samplers/refractive_sampler.jl")
require("data_utils/read_phylosub_data.jl")
require("utils/parallel_utils.jl")
require("utils/general_macros.jl")

#if !isdefined(:filename)
#    println("Using default input file CLL077.csv")    
#    filename = "CLL077.csv"
#end


function run_all_betasplit_experiments(alpha, kind)
    filenames = readdir("../data/phylosub/beta_split")

    alpha_lambda = 1.0
    beta_lambda = 1.0

    for fname in filenames
        if contains(fname, "csv") && contains(fname, kind)
            m = match(r"_([0-9]+).csv", fname)
            index = int(m.captures[1])
            @spawn run_phylo_experiment("beta_split/$fname", alpha, index=index)
        end
    end

end

function run_all_betasplitphylo_experiments(alpha, kind)
    filenames = readdir("../data/phylosub/beta_split_phylo")

    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

    for fname in filenames
        if contains(fname, "csv") && contains(fname, kind)
            job_id, job_ref = smart_spawn(run_phylo_experiment, "beta_split_phylo/$fname", alpha)
        end
    end

end

function run_all_emptysims_experiments(alpha; max_SSMs=Inf)
    filenames = readdir("../data/phylosub/emptysims")

    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

    for fname in filenames
        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", fname)
        N_SSMs = int(m.captures[3])
        if N_SSMs <= max_SSMs
            println("running experiment with $N_SSMs mutations")
            #@spawn run_phylo_experiment("emptysims/$fname", alpha)
           
            job_id, job_ref = smart_spawn( run_phylo_experiment, "emptysims/$fname", alpha)
 

        else
            println("skipping experiment with $N_SSMs mutations")
        end
    end

    

end

function run_aldous_experiments(alpha; max_depth=Inf)
    filenames = readdir("../data/phylosub/aldous")

    for fname in filenames
        m = match(r"\.([0-9]+)\.([0-9]+)\.", fname)
        D = int(m.captures[1])
        if D <= max_depth
            @spawn run_phylo_experiment("aldous/$fname", alpha)
        end
    end
end

function run_phylo_experiment(filename, alpha::Float64;
                              multilocus_filename = nothing, 
                              wl_boundaries::Vector{Float64} = [Inf], #[-1400:50:-1250.0],
                              wl_K_boundaries::Vector{Float64} = [4,5,6,7,Inf],
                              wl_f0::Float64 = 1.0,
                              wl_histogram_test_ratio::Float64 = 0.3,
                              index::Int64 = 0,
                              init_state = nothing )


    trial_index = index

    #model_spec.debug = true

    lambda = 0.2
    gamma = 1.0
    rates_shape = 1.0

    #if filename == "CLL077.csv" || filename == "CLL003.csv"
    #    init_K = 4
    #elseif filename == "CLL006.csv"
    #    init_K = 5
    #end

    if contains(filename, "CLL")
        init_K=4
        filename_base="CLL"
        eta_Temp = 0.0000001
        aisrj_lag = Inf
        num_iterations=10000
    elseif contains(filename, "emptysims")
        filename_base="emptysims"
        eta_Temp = 0.001
        aisrj_lag = Inf
        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", filename)
        init_K = int(m.captures[1])-1
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
        rand_restarts=0
        num_iterations=100000
    elseif contains(filename, "aldous")
        filename_base="aldous"
        eta_Temp = 0.1
        aisrj_lag = Inf
        m = match(r"\.([0-9]+)\.([0-9]+)\.", filename)
        init_K = 4
        D = int(m.captures[1])
        count = int(m.captures[2])

        if D < 100
            wl_K_boundaries -= 1
        end

        rand_restarts=0
        num_iterations=100000
    elseif contains(filename, "betasplit_phylo")
        filename_base="betasplit_phylo"
        aisrj_lag = Inf
        m = match(r"_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)\.", filename)
        init_K = int(m.captures[1])
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
        trial_index = int(m.captures[4])
        rand_restarts=0
        num_iterations=100000
    elseif contains(filename, "betasplit")
        if contains(filename, "chain")
            filename_base = "betasplit_chain"
        elseif contains(filename, "branch")
            filename_base = "betasplit_branch"
        end
        wl_K_boundaries = [3,4,Inf]
        aisrj_lag = Inf
        init_K = 3 
        D = 0
        M_per_cluster = 0
        rand_restarts=0
        num_iterations=100000
    elseif contains(filename, "phylospan")

        if multilocus_filename == nothing
            filename_base = "phylospan"
        else
            filename_base= contains(multilocus_filename, "chain") ? "phylospan_chain" : "phylospan_branch"
        end

        wl_K_boundaries = [3,4,Inf]
        aisrj_lag = Inf
        init_K = 3 
        D = 0
        M_per_cluster = 0
        rand_restarts=0
        num_iterations=10000
    end


    WL_state = WangLandauState(wl_boundaries, wl_K_boundaries, wl_f0, wl_histogram_test_ratio)

    use_latent_rates = true
    #use_latent_rates = false

    @ifndef verbose false

    if !isdefined(:plotting)
        model_spec = ModelSpecification(use_latent_rates, ones(3)/3, false, verbose, false)
    else
        model_spec = ModelSpecification(use_latent_rates, ones(3)/3, false, verbose, plotting)
    end




#    (AA, DD, mu_r, mu_v, names) = read_phylosub_data(filename)
#
#    (M, S) = size(AA)
#    paired_reads = zeros(0,9)
#    if multilocus_filename != nothing
#        name2index = Dict{ASCIIString, Int64}()
#
#        for i = 1:length(names)
#            name2index[names[i]] = i
#        end
#
#        paired_reads = read_multilocus_data(multilocus_filename, name2index)
#    end
#    data = DataState(AA, DD, mu_r, mu_v, paired_reads, names)

    data = constructDataState(filename, multilocus_filename=multilocus_filename)

    result = mcmc(data, lambda, gamma, alpha, rates_shape, init_K, model_spec, num_iterations, 1000, aisrj_lag = aisrj_lag, rand_restarts=rand_restarts, WL_state = WL_state, init_state = init_state)
    (iters, Ks, trainLLs, models) = result

    
    f = open("../results/phylo/$filename_base/$filename_base.ccm.$alpha.$init_K.$D.$M_per_cluster.$trial_index.models", "w")
    serialize(f, models) 
    close(f)
end

