#include("../model/phylo_mcmc.jl")
#include("../data_utils/read_phylosub_data.jl")
#include("../utils/parallel_utils.jl")
#include("../utils/general_macros.jl")

#if !isdefined(:filename)
#    println("Using default input file CLL077.csv")    
#    filename = "CLL077.csv"
#end

phylosub_prefix = "../data/phylosub"
phylospan_prefix = "../data/phylospan"

function run_all_CLL_experiments(filename, alpha)

    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()
    for sample_rate in [0.0001, 0.001, 0.01, 1.0]
        job_id, job_ref = smart_spawn( () -> run_phylo_experiment( filename, alpha, nothing,
                                                                   read_subsample_rate=sample_rate,
                                                                   num_iterations=20000,
                                                                   init_K=4)) 
    end 

end

function run_all_betasplit_experiments(alpha, kind)
    filenames = readdir("../data/phylosub/beta_split")

    alpha_lambda = 1.0
    beta_lambda = 1.0

    for fname in filenames
        if contains(fname, "csv") && contains(fname, kind)
            m = match(r"_([0-9]+).csv", fname)
            index = int(m.captures[1])
            @spawn run_phylo_experiment("$phylosub_prefix/beta_split/$fname", alpha, nothing,
                                         index=index,
                                         num_iterations=100000,
                                         wl_K_boundaries=[3,4,Inf])
        end
    end

end

function run_all_betasplitphylo_experiments(alpha, kind)
    filenames = readdir("../data/phylosub/beta_split_phylo")

    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

    for fname in filenames
        if contains(fname, "csv") && contains(fname, kind)
            
            m = match(r"_([0-9])_", fname)

            init_K = int(m.captures[1])
            job_id, job_ref = smart_spawn( () -> run_phylo_experiment("$phylosub_prefix/beta_split_phylo/$fname", alpha, nothing,
                                                                       num_iterations=100000,
                                                                       init_K=init_K))
        end
    end

end

function run_all_phylospan_experiments(alpha, kind)
    filenames = readdir("../data/phylospan/phylospan_sims")

    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

    base_file = "$phylospan_prefix/phylospan_sims/phylospan.50.ssm.txt"

    for fname in filenames
        if contains(fname, kind)
            job_id, job_ref = smart_spawn( () -> run_phylo_experiment( base_file, alpha, "$phylospan_prefix/phylospan_sims/$fname",
                                                                       num_iterations=10000,
                                                                       wl_K_boundaries=[3,4,Inf] ))
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
          
            init_K = int(m.captures[1])-1
            job_id, job_ref = smart_spawn( () -> run_phylo_experiment( "emptysims/$fname", alpha, nothing,
                                                                        num_iterations=100000,
                                                                        init_K = init_K))
 

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

            wl_K_boundaries = D < 100 ? [3,4,5,6,Inf] : [4,5,6,7,Inf]
            @spawn run_phylo_experiment("aldous/$fname", alpha, nothing,
                                         num_iterations=100000,
                                         wl_K_boundaries=wl_K_boundaries,
                                         init_K=4)
        end
    end
end

# Main entry point

function run_phylo_experiment(filename, alpha::Float64, multilocus_filename; 
                              wl_boundaries::Vector{Float64} = [Inf], #[-1400:50:-1250.0],
                              wl_K_boundaries::Vector{Float64} = [4,5,6,7,Inf],
                              wl_f0::Float64 = 1.0,
                              wl_histogram_test_ratio::Float64 = 0.3,
                              index::Int64 = 0,
                              init_state = nothing,
                              read_subsample_rate = 1.0,
                              num_iterations = 10000,
                              init_K = 3,
                              outputfile=nothing )


    trial_index = index

    #model_spec.debug = true

    lambda = 0.4
    gamma = 1.0
    rates_shape = 1.0

    #if filename == "CLL077.csv" || filename == "CLL003.csv"
    #    init_K = 4
    #elseif filename == "CLL006.csv"
    #    init_K = 5
    #end

    # not useful in conjunction with WL
    rand_restarts=0

    # aisrj too slow
    aisrj_lag = Inf

    multilocus_string = ""

    subsample_rate_string = ""
    if contains(filename, "CLL")
        m = match(r"(CLL[0-9]+)\.csv", filename)
        filename_base=m.captures[1] #"CLL"
        D = 0
        M_per_cluster = 0
        trial_index=1
        subsample_rate_string = ".$read_subsample_rate"
        
    elseif contains(filename, "emptysims")
        filename_base="emptysims"
        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", filename)
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
    elseif contains(filename, "aldous")
        filename_base="aldous"
        m = match(r"\.([0-9]+)\.([0-9]+)\.", filename)
        D = int(m.captures[1])
        count = int(m.captures[2])

    elseif contains(filename, "betasplit_phylo")
        filename_base="betasplit_phylo"
        m = match(r"_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)\.", filename)
        D = int(m.captures[2])
        M_per_cluster = int(m.captures[3])
        trial_index = int(m.captures[4])
    elseif contains(filename, "betasplit")
        if contains(filename, "chain")
            filename_base = "betasplit_chain"
        elseif contains(filename, "branch")
            filename_base = "betasplit_branch"
        end
        D = 0
        M_per_cluster = 0
    elseif contains(filename, "phylospan")

        if multilocus_filename == nothing
            filename_base = "phylospan"
        else
            filename_base= contains(multilocus_filename, "chain") ? "phylospan_chain" : "phylospan_branch"
        end

        m = match(r"\.([0-9]+)\.([0-9]+)\.([0-9]+)\.", multilocus_filename)

        npairs = int(m.captures[1])
        phasing_percent = int(m.captures[2])
        trial_index = int(m.captures[3])
        multilocus_string = "$(npairs)_$phasing_percent"

        D = 0
        M_per_cluster = 0
    end


    WL_state = WangLandauState(wl_boundaries, wl_K_boundaries, wl_f0, wl_histogram_test_ratio)

    #use_latent_rates = true
    use_latent_rates = false

    @ifndef verbose false

    if !isdefined(:plotting)
        model_spec = ModelSpecification(use_latent_rates, ones(3)/3, false, verbose, false)
    else
        model_spec = ModelSpecification(use_latent_rates, ones(3)/3, false, verbose, plotting)
    end



    data = constructDataState(filename, multilocus_filename=multilocus_filename)

    if read_subsample_rate < 1.0
        # subsample in a reproducible way
        r = MersenneTwister(1)
        tied_binomial = (n, n2, p) -> (A = [ rand(r) < p for i = 1:n]; (sum(A), sum(A[1:n2])))

        for i = 1:length(data.reference_counts)
            data.total_counts[i], data.reference_counts[i] = 
                tied_binomial(data.total_counts[i], data.reference_counts[i], read_subsample_rate)
        end
    end

    result = mcmc(data, lambda, gamma, alpha, rates_shape, init_K, model_spec, num_iterations, 1000, aisrj_lag = aisrj_lag, rand_restarts=rand_restarts, WL_state = WL_state, init_state = init_state)
    (iters, Ks, trainLLs, models) = result

    mkpath("../results/phylo/$filename_base")

    if outputfile == nothing
        if contains(filename_base, "phylospan")
            f = open("../results/phylo/$filename_base/$filename_base.$npairs.$phasing_percent.$trial_index.models", "w")
        else
            f = open("../results/phylo/$filename_base/$filename_base.$alpha.$init_K.$D.$M_per_cluster$subsample_rate_string.$trial_index.models", "w")
        end

    else
        try 
            f = open(outputfile, "w")
        catch
            outf = "../results/phylo/$filename_base/$filename_base..$filename..$multilocus_filename.models"
            println("Unable to open $outputfile, outputting to $outf")
            f = open(outf, "w")
        end
    end

    serialize(f, models) 
    close(f)
end

run_phylo_experiment(filename, alpha::Float64) = run_phylo_experiment(filename, alpha, nothing)
