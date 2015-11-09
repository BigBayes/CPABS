#require("phylo_model/phylo_model.jl")
#require("phylo_model/pdf.jl")
#require("phylo_model/phylo_mcmc.jl")
#require("eval/eval_phylo.jl")
#require("data_utils/data_utils.jl")
#require("utils/plot_utils.jl")
#require("utils/probability_util.jl")
#
#using Distributions

function getModelDendrogram(model::ModelState, true_clustering=nothing)
    Z = model.Z
    N::Int = (length(model.tree.nodes) + 1) / 2

    u = zeros(Int64, 2N-1)
    for i=1:length(Z)
        u[Z[i]] += 1
    end 

    ZZ, leaf_times, Etas, inds = model2array(model, return_leaf_times=true)

    if true_clustering != nothing

        A = zeros(Int64, N-1, length(true_clustering))

        for i = 1:length(Z)
            true_cluster = find([any(true_clustering[j] .== i) for j = 1:length(true_clustering)])
           
            A[Z[i]-N, true_cluster] += 1
 
        end

        phis = compute_phis(model)

        annotations = Array(ASCIIString, 2N-1)
        annotations[1:N] = ""
        for i = N+1:2N-1
            annotations[i] = "$(A[i-N,:])\nphi: $(round(phis[i],3))"
        end

        dendrogram(ZZ, u[inds], plot=false, sorted_inds=inds, leaf_times=leaf_times, annotations=annotations)
    else
        dendrogram(ZZ, u[inds], plot=false, sorted_inds=inds, leaf_times=leaf_times)
    end
end

function getModelHierarchy(model::ModelState, true_clustering=nothing, names=nothing)
    Z = model.Z
    N::Int = (length(model.tree.nodes) + 1) / 2

    u = zeros(Int64, 2N-1)
    for i=1:length(Z)
        u[Z[i]] += 1
    end 

    T = GetAdjacencyMatrix(model.tree)

    # really should refactor this
    if true_clustering != nothing || names != nothing

        if true_clustering != nothing
            A = zeros(Int64, N-1, length(true_clustering))
        end
        printed_names = ["" for i = 1:N-1]

        for i = 1:length(Z)
            if true_clustering != nothing
                true_cluster = find([any(true_clustering[j] .== i) for j = 1:length(true_clustering)])
                A[Z[i]-N, true_cluster] += 1
            end

            if names != nothing 
                old_name = printed_names[Z[i]-N] 
                printed_names[Z[i]-N] = old_name == "" ? names[i] : "$old_name, $(names[i])"
            end
 
        end

        phis = compute_phis(model)
        perm = reverse(sortperm(sum(phis[N+1:end,:],2)[:]))
        perm_phis = phis[N+1:end,:][perm,:]

        annotations = Dict{Int64, ASCIIString}()
        for i = N+1:2N-1
            if names != nothing
                annotations[i-N] = "$(printed_names[perm[i-N]])\nphi: $([round(perm_phis[i-N,j],3) for j = 1:length(perm_phis[i-N,:])])"
            else
                annotations[i-N] = "$(A[perm[i-N],:])   phi: $(round(mean(perm_phis[i-N,:]),3))"
            end
        end

        plot_subclonal_hierarchy(T[perm,perm], annotations)
    else
        plot_subclonal_hierarchy(T, nothing)
    end
end

function getModels(models_filename::ASCIIString)
    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)
    models
end

function genModelsNavigator(models_filename::ASCIIString; true_clustering=nothing)

    true_clustering = true_clustering
    models = getModels(models_filename)

    curindex = 1 

    function showModel(index::Int64)
        p = getModelHierarchy(models[index], true_clustering)
        curindex = index
        return p
    end

    function stepN(N::Int64)
        curindex += N
        curindex = curindex > length(models) ? length(models) : 
                   curindex < 1 ? 1 : curindex
        showModel(curindex) 
    end 

    showModel, stepN
end

#function computeLogPosteriors(models, data::DataState)
#    latent_rates = models[1].rates != zeros(length(models[1].rates))
#    model_spec = ModelSpecification(latent_rates, zeros(3), false, false, false)
#
#    wl_state = models[end].WL_state
#    total_partition_mass = logsumexp(wl_state.partition_function)
#
#    log_posteriors = zeros(length(models))
#
#    for n = 1:length(models)
#        model = models[n]
#        N::Int = (length(model.tree.nodes)+1)/2
#
#        log_p = full_pdf(model, model_spec, data)
#        w = get_partition_function( wl_state, N-1, log_p)
#        w = exp(w - total_partition_mass)
#
#
#
#    end
#
#    
#end

function getRepresentativeSample(models; sampled_bin=nothing)
  
    wl_state = models[end].WL_state 
    @assert length(wl_state.energy_boundaries) == 1

    Ns = [int((length(models[i].tree.nodes)+1)/2) for i = 1:length(models)]

    bin_map = [get_bin(wl_state, N-1, 0.0)[1] for N = 1:maximum(Ns)]

    bins = [bin_map[N] for N in Ns] 

    bin_logprobs = wl_state.partition_function[:,1]
    sampled_bin = sampled_bin == nothing ? rand(Categorical(exp_normalize(bin_logprobs))) : sampled_bin

    I = find(bins .== sampled_bin)
    model_index = I[rand(1:length(I))]

    models[model_index]
end

function plotRepresentativeSample(models; true_clustering=nothing, bin=nothing, data=nothing, names=nothing)
    model = getRepresentativeSample(models, sampled_bin=bin)
    p = getModelHierarchy(model, true_clustering, names)

    if data != nothing
        model_spec = ModelSpecification(false, zeros(3), false, false, false)
        LL = likelihood(model, model_spec, data)
        PR = prior(model, model_spec)
        println("Model prior, log-likelihood: $PR, $LL")
    end

    p
end 

function getWLPartitionFuncs(models)

    K = length(models[1].WL_state.partition_function)
    N = length(models)
   
    Zs = zeros(K,N) 

    for i = 1:N

        model = models[i]
        Z = model.WL_state.partition_function

        Zs[:,i] = Z
    end

    Zs
end

function plotModelClusteringMatrices(models_filename::ASCIIString)
    models = getModels(models_filename)

    fname = models_filename

    m = match(r".*/([^/]*).ccm\.([0-9]+)\..*\.([0-9]+[0-9]*)\.models", fname)
    fname_base = m.captures[1] 
#    alpha = float(m.captures[2]) 
    index = int(m.captures[3]) 
    data_file = "beta_split/$(fname_base)_$index.csv"
    clusters_file = "$(fname_base)_clusters_$index.jld"

    data = constructDataState(data_file)
    clusters = load("../data/phylosub/beta_split/$clusters_file", "clusters")

    predicted_clustering = compute_cocluster_matrix(models_filename, data)
    predicted_ancestorship = compute_ancestorship_matrix(models_filename, data)

    t = Table(2,2)

    hold(true)

    p1 = FramedPlot(title="Predicted Ancestorship")
    p2 = FramedPlot(title="Predicted Coclustering")

    Winston.ghf(p1)
    imagesc(predicted_ancestorship)
    Winston.ghf(p2)
    imagesc(predicted_clustering)

    t[1,1] = p1
    t[1,2] = p2

    p3 = FramedPlot(title="True Ancestorship")
    p4 = FramedPlot(title="True Coclustering")

    ground_truth_ancestorship = get_ancestorship_from_clusters(clusters, "chain")
    ground_truth_clustering = get_coclustering_from_clusters(clusters)

    Winston.ghf(p3)
    imagesc(ground_truth_ancestorship)
    Winston.ghf(p4)
    imagesc(ground_truth_clustering)

    t[2,1] = p3
    t[2,2] = p4

    display(t)
    t
end

function plot_phi_v_time(N, M, alpha)

    X = Float64[]
    Y = Float64[] 

    for m = 1:M
        model = draw_random_tree(N, N, 1, 1.0, 1.0, alpha, 1.0, WangLandauState()) 

        tree = model.tree

        times = compute_times(model)
        phis = compute_phis(model)

        for i = N+1:2N-1
            push!(X, times[i])
            push!(Y, phis[i])
        end
    end

    pts = Points(X,Y, kind="dot")
    p = FramedPlot(title="alpha= $alpha", xlabel="time", ylabel="phi")
    add(p, pts)
    X,Y,p
end 
