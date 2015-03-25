require("phylo_model/phylo_model.jl")
require("eval/eval_phylo.jl")
require("utils/plot_utils.jl")

function getModelDendrogram(model::ModelState)
    Z = model.Z
    N::Int = (length(model.tree.nodes) + 1) / 2

    u = zeros(Int64, 2N-1)
    for i=1:length(Z)
        u[Z[i]] += 1
    end 

    ZZ, leaf_times, Etas, inds = model2array(model, return_leaf_times=true)

    dendrogram(ZZ, u[inds], plot=false, sorted_inds=inds, leaf_times=leaf_times)
end

function genModelsNavigator(models_filename::ASCIIString)
    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    curindex = 1 

    function showModel(index::Int64)
        p = getModelDendrogram(models[index])
        display(p)
        curindex = index
    end

    function stepN(N::Int64)
        curindex += N
        curindex = curindex > length(models) ? length(models) : 
                   curindex < 1 ? 1 : curindex
        showModel(curindex) 
    end 

    showModel, stepN
end

function plotModelClusteringMatrices(models_filename::ASCIIString)
    f = open(models_filename, "r")
    models = deserializeModels(f)
    close(f)

    fname = models_filename

    m = match(r"(.*).ccm\.([0-9]+)\..*\.([0-9]+[0-9]*)\.models", fname)
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

    p1 = FramedPlot()
    p2 = FramedPlot()

    Winston.ghf(p1)
    imagesc(predicted_ancestorship)
    Winston.ghf(p2)
    imagesc(predicted_clustering)

    t[1,1] = p1
    t[1,2] = p2

    p3 = FramedPlot()
    p4 = FramedPlot()

    ground_truth_ancestorship = get_ancestorship_from_clusters(clusters, "chain")
    ground_truth_clustering = get_coclustering_from_clusters(clusters)

    Winston.ghf(p3)
    imagesc(ground_truth_ancestorship)
    Winston.ghf(p4)
    imagesc(ground_truth_clustering)

    t[2,1] = p3
    t[2,2] = p4

    display(t)
end
