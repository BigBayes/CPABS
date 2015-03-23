require("phylo_model/phylo_model.jl")
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
