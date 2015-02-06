using HDF5, JLD

require("utils/plot_utils.jl")
require("eval/eval_utils.jl")
require("phylo_model/phylo_model.jl")
require("data_utils/read_phylosub_data.jl")
require("phylo_model/phylo_mcmc.jl")

function eval_phylo_experiments(path, filename_base; p=nothing, offset= 0.0, color="red", label="", kind="diamond", alpha=1)
    filenames = readdir(path)

    delim = contains(filename_base, "pwgs") ? ' ' : ','

    if contains(filename_base, "aldous")
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
        
        S = Array(Any, length(n_clusters), length(depths), length(n_mutations))

        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".models")
                m = match(r"\.1\.0\.([0-9]+)\.([0-9]+)\.([0-9]+)\.models", fname)
                C = int(m.captures[1])
                D = int(m.captures[2])
                N = int(m.captures[3])

                

                if N > 500
                    continue
                end

                data_file = "emptysims/emptysim.$C.$D.$N.ssm.txt"
                Ytrue = get_true_clustering_emptysims(C, N)
                models_fname = "$path/$fname"

                

                C_index = find(C .== n_clusters)[1]
                D_index = find(D .== depths)[1]
                N_index = find(N .== n_mutations)[1]

                println("Beginning evaluation ($C, $D, $N) from .models")
                S[C_index, D_index, N_index] = @spawn eval_emptysims_experiment(data_file, models_fname, Ytrue)

            elseif contains(fname, filename_base) && contains(fname, ".csv")
                m = match(r"\.1\.0\.([0-9]+)\.([0-9]+)\.([0-9]+)\.csv", fname)
                C = int(m.captures[1])
                D = int(m.captures[2])
                N = int(m.captures[3])

                m = match(r"(.*)\.csv", fname)
                fname_base = m.captures[1]
      
                # skip if we have the new .models runs 
                if N > 500 || "$fname_base.models" in filenames
                    continue
                end

                triu_inds = find(triu(ones(C*N,C*N),1))

                C_index = find(C .== n_clusters)[1]
                D_index = find(D .== depths)[1]
                N_index = find(N .== n_mutations)[1]

                Ypred = readdlm("$path/$fname", delim)
                Ytrue = get_true_clustering_emptysims(C, N)

                println("Beginning evaluation ($C, $D, $N) from .csv")
                S[C_index, D_index, N_index] = @spawn aupr(Ypred[triu_inds], Ytrue[triu_inds])
            end
        end

        println("Fetching results")
        for i = 1:length(S)
            auprs[i] = fetch(S[i])
        end

        return auprs
    elseif contains(filename_base, "betasplit")

        auprs = zeros(8)
        cauprs = zeros(8)
        aucs = zeros(8)

        S = Array(Any,8)
        for fname in filenames
            if contains(fname, filename_base) && contains(fname, ".models")
                m = match(r"(.*).ccm\.([0-9]+)\..*\.([0-9]+)\.models", fname)
                fname_base = m.captures[1] 
                alpha = float(m.captures[2]) 
                index = int(m.captures[3]) 
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
        for i = 1:8
            auprs[i], cauprs[i] = fetch(S[i]) #eval_betasplit_experiment(data_file, models_fname, clusters, kind) 
        end

        return auprs, cauprs
    end
end



function eval_phylospan_experiment(data_filename, multilocus_filename, models_filename, kind)

    data = constructDataState(data_filename, multilocus_filename=multilocus_filename)

    predicted = compute_ancestorship_matrix(models_filename, data)
    ground_truth = get_true_ancestorship_phylospan(kind)
  
    N, N = size(ground_truth)
 
    nondiagonals = find(1 .- eye(N)) 

    aupr(predicted[nondiagonals], ground_truth[nondiagonals])
end

function eval_betasplit_experiment(data_filename, models_filename, clusters, kind)

    data = constructDataState(data_filename)

    predicted = compute_ancestorship_matrix(models_filename, data)
    ground_truth = get_ancestorship_from_clusters(clusters, kind)
  
    N, N = size(ground_truth)
 
    nondiagonals = find(1 .- eye(N)) 

    a_aupr = aupr(predicted[nondiagonals], ground_truth[nondiagonals])


    predicted = compute_cocluster_matrix(models_filename, data)
    ground_truth = get_coclustering_from_clusters(clusters)

    triu_inds = find(triu(ones(N,N),1))
    c_aupr = aupr(predicted[triu_inds], ground_truth[triu_inds])

    a_aupr, c_aupr
end

function eval_emptysims_experiment(data_filename, models_filename, ground_truth)

    data = constructDataState(data_filename)

    predicted = compute_cocluster_matrix(models_filename, data)

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
    
    return Y
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

function compute_cocluster_matrix(models_filename::ASCIIString, data::DataState)

    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    wl_state = models[end].WL_state
    M = length(models[end].Z)
    
    model_spec = ModelSpecification(zeros(3), false, false, false)

    total_partition_mass = logsumexp(wl_state.partition_function)

    cocluster_matrix = zeros(M,M)
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
    

    return cocluster_matrix
end

function compute_ancestorship_matrix(models_filename::ASCIIString, data::DataState)

    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    wl_state = models[end].WL_state
    M = length(models[end].Z)

    total_partition_mass = logsumexp(wl_state.partition_function)

    model_spec = ModelSpecification(zeros(3), false, false, false)

    ancestorship_matrix = zeros(M,M)
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
        aucs[C_index, D_index, N_index] = float(A[i,5])
    end

    auprs, aucs
end
