using HDF5, JLD

require("eval/eval_utils.jl")
require("phylo_model/phylo_model.jl")
require("data_utils/read_phylosub_data.jl")
require("phylo_model/phylo_mcmc.jl")
require("utils/parallel_utils.jl")

function eval_phylo_experiments(path, filename_base; p=nothing, offset= 0.0, color="red", label="", kind="diamond", alpha=1)
    filenames = readdir(path)

    delim = contains(filename_base, "pwgs") ? ' ' : ','

    if contains(filename_base, "CLL")
        m = match(r"CLL([0-9]+)", filename_base)
        CLL_id = m.captures[1] 
        clusters, Ytrue= get_true_clustering_CLL(CLL_id)

        subsample_rates = [0.0001, 0.001, 0.01, 1.0]

        auprs = zeros(length(subsample_rates))
        S = Array(Any, length(subsample_rates))

        smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

        for fname in filenames
            if contains(fname, filename_base) && contains(fname, "CLL$CLL_id")
                m = match(r"\.([01]\.[0-9]+)\.[0-9]\.models", fname)
                subsample_rate = float(m.captures[1])

                rate_index = find(subsample_rates .== subsample_rate)[1]
                data_file = "../data/phylosub/CLL$CLL_id.csv"
                models_fname = "$path/$fname"

                job_id, S[rate_index] = smart_spawn( () -> eval_CLL_experiment(data_file, models_fname, Ytrue, subsample_rate))

            end
        end

        println("Fetching results")
        for i = 1:length(S)
            auprs[i] = fetch(S[i])
        end

        return auprs

    elseif contains(filename_base, "aldous")
        Ytrue = get_true_clustering()
        depths = [20,30,50,70,100,200,300]
        auprs = zeros(length(depths), 10)
        aucs = zeros(length(depths), 10)

        triu_inds = find(triu(ones(288,288),1))


        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".csv")
                m = match(r"\.([0-9]+)\.([0-9]+)\.csv", fname)
                D = int(m.captures[1])
                count = int(m.captures[2])

                D_index = find(D .== depths)[1]

                Ypred = readdlm("$path/$fname", delim)
                auprs[D_index, count] = aupr(Ypred[triu_inds], Ytrue[triu_inds])
            end
        end
        meds = median(auprs,2)
        mins = minimum(auprs,2)
        maxs = maximum(auprs,2)

        if p == nothing
            p = FramedPlot()
        end

        ebars = Winston.ErrorBarsY(depths+offset, mins, maxs, color=color)
        add(p, ebars)
        pts = Points(depths+offset, meds, color=color, size=0.5, kind=kind)
        setattr(pts, label=label)
        add(p, pts)
        display(p)

        return auprs, aucs, p, pts

    elseif contains(filename_base, "emptysim")

        n_clusters = [2,3,4,5]
        depths = [20,30,50,70,100,200,300]
        n_mutations = [5,10,25,50,100,200,500]


        auprs = zeros(length(n_clusters), length(depths), length(n_mutations))
        
        S = Array(RemoteRef, length(n_clusters), length(depths), length(n_mutations))

        spawn_counter = 0
        running_procs = Array(Any, nworkers())

        available_workers = workers()

        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".models")
		println("evaluating $fname")
                m = match(r"\.1\.0\.([0-9]+)\.([0-9]+)\.([0-9]+)\.([0-9]*).?models", fname)
		@assert m != Nothing
                C = int(m.captures[1])
                D = int(m.captures[2])
                N = int(m.captures[3])

                

                if N > 500
                    continue
                end

                data_file = "emptysims/emptysim.$(C+1).$D.$N.ssm.txt"
                Ytrue = get_true_clustering_emptysims(C, N)
                models_fname = "$path/$fname"

                

                C_index = find(C .== n_clusters)[1]
                D_index = find(D .== depths)[1]
                N_index = find(N .== n_mutations)[1]

                

                s = RemoteRef()
                worker = pop!(available_workers)
                println("Beginning evaluation ($C, $D, $N) on worker $worker")
                @async put!(s, remotecall_fetch( worker, eval_emptysims_experiment, data_file, models_fname, Ytrue) )
                S[C_index, D_index, N_index] = s
                running_procs[worker-1] = S[C_index, D_index, N_index]

            end

            while length(available_workers) == 0
                sleep(100)
                for i = 1:nworkers()
                    if isready(running_procs[i])
                        push!(available_workers, i+1)
                    end
                end
            end
        end

        println("Fetching results")
        for i = 1:length(S)
            auprs[i] = fetch(S[i])
        end

        return auprs
    elseif contains(filename_base, "betasplit_phylo")

        smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

        n_clusters = [3,4,5]
        depths = [50,70,100]
        n_mutations = [10, 25, 100]

        auprs = zeros(3,3,3,8)
        cauprs = zeros(3,3,3,8)


        precs = Array(Any, 3,3,3,8)
        recalls = Array(Any, 3,3,3,8)
        c_precs = Array(Any, 3,3,3,8)
        c_recalls = Array(Any, 3,3,3,8)


        S = Array(Any,3,3,3,8)
        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".models") && !contains(fname, "ccm")
                m = match(r"(.*)(\.ccm)?\.0\.0\.([0-9]+)\.([0-9]+[0-9]*)\.([0-9]+[0-9]*)\.([0-9]+)\.models", fname)
                fname_base = m.captures[1] 

                C = int(m.captures[3])
                D = int(m.captures[4])
                N = int(m.captures[5])
                index = int(m.captures[6])

                C_index = find(C .== n_clusters)[1]
                D_index = find(D .== depths)[1]
                N_index = find(N .== n_mutations)[1]

                data_file = "../data/phylosub/beta_split_phylo/$(fname_base)_chain_$(C)_$(D)_$(N)_$index.csv"
                clusters_file = "gt_$(fname_base)_chain_$(C)_$(D)_$(N)_$index.jld"

                clusters = load("../data/phylosub/beta_split_phylo/$clusters_file", "clusters")

                models_fname = "$path/$fname"
                kind = "chain" 

                println("evaluating trial $C $D $N $index")
                #S[index] = @spawn eval_betasplit_experiment(data_file, models_fname, clusters, kind) 
                job_id, S[C_index, D_index, N_index, index] = smart_spawn(eval_betasplit_experiment, data_file, models_fname, clusters, kind) 
            end

        end
        for i = 1:length(S)
            auprs[i], precs[i], recalls[i], cauprs[i], c_precs[i], c_recalls[i] = fetch(S[i])
        end

        return auprs, precs, recalls, cauprs, c_precs, c_recalls
    elseif contains(filename_base, "betasplit")

        auprs = zeros(40)
        cauprs = zeros(40)
        aucs = zeros(40)

        S = Array(Any,40)
        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".models")
                m = match(r"(.*)(\.ccm)?\.([0-9]+)\..*\.([0-9]+[0-9]*)\.models", fname)
                fname_base = m.captures[1] 
                alpha = float(m.captures[3]) 
                index = int(m.captures[4]) 
                data_file = "beta_split/$(fname_base)_$index.csv"
                clusters_file = "$(fname_base)_clusters_$index.jld"

                clusters = load("../data/phylosub/beta_split/$clusters_file", "clusters")

                models_fname = "$path/$fname"
                kind = contains(fname, "chain") ? "chain" : 
                       contains(fname, "branch") ? "branch" : nothing

                println("evaluating trial $index")
                S[index] = @spawn eval_betasplit_experiment(data_file, models_fname, clusters, kind) 
            end

        end
        for i = 1:40
            auprs[i], cauprs[i] = fetch(S[i]) #eval_betasplit_experiment(data_file, models_fname, clusters, kind) 
        end

        return auprs, cauprs
    end
end

function eval_phylospan_experiments(path, filename_base)
    filenames = readdir(path)
    
    smart_spawn, get_job_counter, get_jobs = initialize_smart_spawn()

    data_filename = "../data/phylospan/phylospan_sims/phylospan.50.ssm.txt"

    npairs = [5,10,20,50,100]
    phasing_percents = [0,10,20,30,40,50]
    
    auprs = zeros(5,6,10)

    S = Array(Any, 5,6,10)

    for fname in filenames
        if contains(fname, filename_base)
            m = match(r".*\.([0-9]+)\.([0-9]+)\.([0-9]+)\.models", fname)

            npair = int(m.captures[1])
            phasing_percent = int(m.captures[2])
            index = int(m.captures[3])

            multilocus_filename = "../data/phylospan/phylospan_sims/phylospan.$filename_base.$npair.$phasing_percent.$index.ssm.txt"

            
            npair_index = find(npairs .== npair)[1]
            phasing_index = find(phasing_percents .== phasing_percent)[1]

            println("evaluating experiment $npair_index, $phasing_index, $index")  

            job_id, S[npair_index, phasing_index, index] = smart_spawn(eval_phylospan_experiment, data_filename, multilocus_filename, "$path/$fname", filename_base) 

        end
    end
    
    for i = 1:length(S)
        auprs[i] = fetch(S[i])
    end

    return auprs
end

function eval_phylospan_experiment(data_filename, multilocus_filename, models_filename, kind)

    data = constructDataState(data_filename, multilocus_filename=multilocus_filename)

    predicted = compute_ancestorship_matrix(models_filename, data)
    ground_truth, ground_truth_cluster = get_true_ancestorship_phylospan(kind)
  
    N, N = size(ground_truth)

    triu_inds = find(triu(ones(N,N),1))
    tril_inds = find(tril(ones(N,N),-1))
 
    N_a = sum(ground_truth[triu_inds])
    N_d = sum(ground_truth[tril_inds])

    a_aupr = aupr(predicted[triu_inds], ground_truth[triu_inds])
    d_aupr = aupr(predicted[tril_inds], ground_truth[tril_inds])
   
    predicted_cluster = compute_cocluster_matrix(models_filename, data)

    N_c = sum(ground_truth_cluster[triu_inds])
    c_aupr = aupr(predicted_cluster[triu_inds], ground_truth_cluster[triu_inds])

    P_truth = ground_truth_cluster + ground_truth + ground_truth'
    @assert length(find(P_truth[triu_inds] .> 1)) == 0
    P_truth[find(P_truth .!=0)] = 1.0
    P_truth = 1.0 - P_truth

    P_pred = predicted_cluster + predicted + predicted'
    P_pred = 1.0 - P_pred

    n_aupr = aupr(P_pred[triu_inds], P_truth[triu_inds])

    N_n = length(find(P_truth[triu_inds] .== 1))


    (a_aupr * N_a + d_aupr * N_d + c_aupr * N_c + n_aupr * N_n) / (N_a + N_d + N_c + N_n)
end

function eval_betasplit_experiment(data_filename, models_filename, clusters, kind)

    data = constructDataState(data_filename)

    predicted = compute_ancestorship_matrix(models_filename, data, start_index=50000)
    ground_truth = get_ancestorship_from_clusters(clusters, kind)
  
    N, N = size(ground_truth)
 
    nondiagonals = find(1 .- eye(N)) 

    a_aupr, a_precs, a_recalls = aupr(predicted[nondiagonals], ground_truth[nondiagonals], return_curves=true)


    predicted = compute_cocluster_matrix(models_filename, data, start_index=50000)
    ground_truth = get_coclustering_from_clusters(clusters)

    triu_inds = find(triu(ones(N,N),1))
    c_aupr, c_precs, c_recalls = aupr(predicted[triu_inds], ground_truth[triu_inds], return_curves=true)

    a_aupr, a_precs, a_recalls, c_aupr, c_precs, c_recalls
end

function eval_emptysims_experiment(data_filename, models_filename, ground_truth)

    data = constructDataState(data_filename)

    predicted = compute_cocluster_matrix(models_filename, data)
    N, N = size(predicted)

    triu_inds = find(triu(ones(N,N),1))
    c_aupr = aupr(predicted[triu_inds], ground_truth[triu_inds])

    c_aupr
end

function eval_CLL_experiment(data_filename, models_filename, ground_truth, subsample_rate)

    data = constructDataState(data_filename)
    if subsample_rate < 1.0
        # subsample in a reproducible way
        r = MersenneTwister(1)
        tied_binomial = (n, n2, p) -> (A = [rand(r) < p for i = 1:n]; (sum(A), sum(A[1:n2])))

        for i = 1:length(data.reference_counts)
            data.total_counts[i], data.reference_counts[i] = 
                tied_binomial(data.total_counts[i], data.reference_counts[i], subsample_rate)
        end
    end

    predicted = compute_cocluster_matrix(models_filename, data)
    N, N = size(predicted)

    triu_inds = find(triu(ones(N,N),1))
    c_aupr = aupr(predicted[triu_inds], ground_truth[triu_inds])

    c_aupr

end

function get_true_clustering()
    Y = zeros(288,288)
    ranges = { 1:108, 109:192, 193:252, 253:288 }
    for k = 1:4
        Y[ranges[k], ranges[k]] = 1
    end
    Y
end

function get_true_clustering_emptysims(C, N)
    base_range = 1:N
    ranges = { base_range + i*N for i = 0:C-1}

    Y = zeros(C*N,C*N)
    for k = 1:C
        Y[ranges[k], ranges[k]] = 1
    end
    Y
end

function get_true_ancestorship_phylospan(kind)

    Y = zeros(300,300)

    if kind == "chain"
        Y[1:100, 101:end] = 1
        Y[101:200, 201:end] = 1    
    elseif kind == "branch"
        Y[1:100, 101:end] = 1
    else
        Y[1:100, 101:end] = 1
        Y[101:200, 201:end] = 0.5
    end
   
    Y2 = zeros(300,300)
   
    Y2[1:100, 1:100] = 1 
    Y2[101:200, 101:200] = 1 
    Y2[201:300, 201:300] = 1 
 
    return Y, Y2
end

#      CLL077      CLL003      CLL006
# 1    BCL2L13     ADAD1       ARHGAP29
# 2    COL24A1     AMTN        EGFR
# 3    DAZAP1      APBB2       IRF4
# 4    EXOC6B      ASXL1       KIAA0182
# 5    GHDC        ATM         KIAA0319L
# 6    GPR158      BPIL2       KLHL4
# 7    HMCN1       CHRNB2      MED12
# 8    KLHDC2      CHTF8       PILRB
# 9    LRRC16A     FAT3        RBPJ
# 10   MAP2K1      HERC2       SIK1
# 11   NAMPTL      IL11RA      U2AF1
# 12   NOD1        MTUS1
# 13   OCA2        MUSK
# 14   PLA2G16     NPY
# 15   SAMHD1      NRG3
# 16   SLC12A1     PLEKHG5
# 17               SEMA3E
# 18               SF3B1
# 19               SHROOM1
# 20               SPTAN1



function get_true_clustering_CLL(id)
    if id == "077"
        groups = { [1,11,6,15,16], [3,4,5,13,14], [9], [8,2,12,7,10] } 
    elseif id == "003"
        groups = { [1,8,10,11,18,19], [4,13,17], [14,15,9,7], [2,3,5,6,16,12,20] }
    elseif id == "006"
        groups = { [1,2,4,6,7,8,10], [11], [5], [3], [9] }
    end

    N = maximum([maximum(g) for g in groups])
    Y = zeros(N,N)

    for g in groups
        Y[g, g] = 1
    end

    return (groups, Y)
end

# assumes clusters is ordered from root to leaves
function get_ancestorship_from_clusters(clusters, kind)
    Ns = [length(clusters[i]) for i = 1:length(clusters)]
    N = sum(Ns)

    Y = zeros(N,N)

    if kind == "chain"
        for i = 1:length(clusters)
            for j = i+1:length(clusters)
                Y[clusters[i], clusters[j]] = 1
            end
        end        
    elseif kind == "branch"
        for i = 2:length(clusters)
            Y[clusters[1], clusters[i]] = 1
        end
    else
        error("not implemented")
    end

    return Y
end

function get_coclustering_from_clusters(clusters)
    Ns = [length(clusters[i]) for i = 1:length(clusters)]
    N = sum(Ns)

    Y = zeros(N,N)
   
    for i = 1:length(clusters)
        Y[clusters[i], clusters[i]] = 1
    end 

    return Y
end

function compute_cocluster_aupr(predicted, ground_truth; return_curves=false)
    N, N = size(ground_truth)
 
    nondiagonals = find(1 .- eye(N)) 
    aupr(predicted[nondiagonals], ground_truth[nondiagonals], return_curves=return_curves)
end

function plot_paired_scatter(auprs1, auprs2, label1, label2, row_index, title)

    color1 = "red"
    color2 = "blue"

    kind1 = "filled circle"
    kind2 = "filled square"

    p = FramedPlot()

    c1 = Points([1:10], auprs1[row_index, :], kind=kind1, color=color1)
    c2 = Points([1:10], auprs2[row_index, :], kind=kind2, color=color2)
   
    setattr(c1, label=label1)
    setattr(c2, label=label2)

    add(p,c1)
    add(p,c2)

    min_aupr = minimum([auprs1[row_index,:] , auprs2[row_index, :]])

    max_aupr = maximum([auprs1[row_index,:] , auprs2[row_index, :]])
    
    setattr(p, yrange=( min_aupr - 0.15, max_aupr+0.03))

    l = Legend(0.65, 0.2, {c1, c2})

    add(p, l)

    setattr(p, xlabel="Trial")
    setattr(p, ylabel="AUPR")
    setattr(p, title=title)
    display(p)

    p    
end

function compute_cocluster_matrix(models_filename::ASCIIString, data::DataState; start_index=1, latent_rates=false, K=nothing)

    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    wl_state = models[end].WL_state
    M = length(models[end].Z)
   
    #latent_rates = models[1].rates != zeros(length(models[1].rates))
 
    model_spec = ModelSpecification(latent_rates, zeros(3), false, false, false)

    total_partition_mass = logsumexp(wl_state.partition_function)

    cocluster_matrix = zeros(M,M)
    sum_w = 0.0
    for n = start_index:length(models)
        model = models[n]
        N::Int = (length(model.tree.nodes)+1)/2

        if K != nothing && N-1 != K
            continue
        end

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
    

    return cocluster_matrix
end

function compute_ancestorship_matrix(models_filename::ASCIIString, data::DataState; start_index=1, latent_rates=false)

    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    wl_state = models[end].WL_state
    M = length(models[end].Z)

    total_partition_mass = logsumexp(wl_state.partition_function)

    #latent_rates = models[1].rates != zeros(length(models[1].rates))
    model_spec = ModelSpecification(latent_rates, zeros(3), false, false, false)

    ancestorship_matrix = zeros(M,M)
    sum_w = 0.0
    for n = start_index:length(models)
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
        for c1_index in 1:length(C)
            for c2_index in 1:length(C)
                if IsRightAncestor(model.tree, c1_index, c2_index) && c1_index != c2_index
                    c1 = C[c1_index]
                    c2 = C[c2_index]
                    ancestorship_matrix[c1, c2] += w
                end 
            end
        end
    end
    ancestorship_matrix /= sum_w
    

    return ancestorship_matrix
end

function read_phylosub_results()
    A, header = readdlm("../results/emsim.results.txt",'\t', header=true)

    n_clusters = [2,3,4,5]
    depths = [20,30,50,70,100,200,300]
    n_mutations = [5,10,25,50,100,200,500]

    auprs = zeros(length(n_clusters), length(depths), length(n_mutations))

    for i = 1:size(A,1)
        C = int(A[i,1])-1
        D = int(A[i,2])
        N = int(A[i,3])

        if N > 500
            continue
        end

        C_index = find(C .== n_clusters)[1]
        D_index = find(D .== depths)[1]
        N_index = find(N .== n_mutations)[1]

        auprs[C_index, D_index, N_index] = float(A[i,4])
        #aucs[C_index, D_index, N_index] = float(A[i,5])
    end

    auprs #, aucs
end

function eval_betasplit_phylo_results()
    path = "../results/phylo/betasplit_phylo"

    n_clusters = [3,4,5]
    depths = [50,70,100]
    n_mutations = [10, 25, 100]

    auprs = zeros(3,3,3,8)
    precs = Array(Any, 3,3,3,8)
    recalls = Array(Any, 3,3,3,8)

    a_auprs = zeros(3,3,3,8)
    a_precs = Array(Any, 3,3,3,8)
    a_recalls = Array(Any, 3,3,3,8)
 
    filenames = readdir(path)

    for fname in filenames
        m = match(r"betasplit_phylo_chain_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+).csv.*", fname)
        fname_base = "betasplit_phylo"

        C = int(m.captures[1])
        D = int(m.captures[2])
        N = int(m.captures[3])
        index = int(m.captures[4])

        C_index = find(C .== n_clusters)[1]
        D_index = find(D .== depths)[1]
        N_index = find(N .== n_mutations)[1]

        clusters_file = "gt_$(fname_base)_chain_$(C)_$(D)_$(N)_$index.jld"
        clusters = load("../data/phylosub/beta_split_phylo/$clusters_file", "clusters")
        ground_truth = get_coclustering_from_clusters(clusters)

        phylosub_results = readdlm("$path/$fname", ' ')

        predicted_lt = phylosub_results[1,:][:]

        predicted = zeros(C*N,C*N)
        predicted[find(tril(ones(C*N,C*N),-1))] = predicted_lt
        predicted = predicted + predicted'
       
        area, precision, recall = compute_cocluster_aupr(predicted, ground_truth, return_curves=true)

        auprs[C_index, D_index, N_index, index] = area
        precs[C_index, D_index, N_index, index] = precision
        recalls[C_index, D_index, N_index, index] = recall 

        ancestor_lt = phylosub_results[2,:][:]
        descendant_lt = phylosub_results[3,:][:]

        predicted = zeros(C*N,C*N)
        predicted[find(tril(ones(C*N,C*N),-1))] = descendant_lt
        predicted = predicted'
        predicted[find(tril(ones(C*N,C*N),-1))] = ancestor_lt

        ground_truth = get_ancestorship_from_clusters(clusters, "chain")

        nondiagonals = find(1 .- eye(C*N)) 

        area, precision, recall = aupr(predicted[nondiagonals], ground_truth[nondiagonals], return_curves=true)
 
        a_auprs[C_index, D_index, N_index, index] = area
        a_precs[C_index, D_index, N_index, index] = precision
        a_recalls[C_index, D_index, N_index, index] = recall 
    end

    a_auprs, a_precs, a_recalls, auprs, precs, recalls 
end
