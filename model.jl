load("tree.jl")

type ModelState
    lambda::Float64
    gamma::Float64
    w_sigma::Float64
    tree::Tree{Int64}
    weights::Array{Float64,2}
    beta::Array{Float64,1}
    beta_p::Array{Float64,1}
    beta_c::Array{Float64,1}
    a::Array{Float64,1}
    b::Array{Float64,1}
    c::Float64
end

# Tree/Model Utility Functions

# Prune tree while also adjusting the weight matrix appropriately
function prune_tree!(model::ModelState,
                     prune_index::Int)
    tree = model.tree
    weight_indices = weight_index_pointers(tree)

    parent = tree.nodes[prune_index].parent

    if parent.children[0].index == prune_index
        sibling = parent.children[1]
    else
        sibling = parent.children[0]
    end


    sibling_start = weight_indices[sibling.index]
    parent_start = weight_indices[parent.index]

    sibling_end = sibling_start + sibling.location - 1
    parent_end = parent_start + parent.location - 1

    weight_permutation = linspace(1,size(model.weights)[1],size(model.weights)[1])


    if sibling_start < parent_start
        del(weight_permutation, parent_start:parent_end)
    
        for p = reverse(parent_start:parent_end)
            insert(weight_permutation, start_end + 1, p)
        end
    else
        for p = reverse(parent_start:parent_end)
            insert(weight_permutation, start_end + 1, p)
        end

        del(weight_permutation, parent_start:parent_end)
    end

    permute_rows_and_cols(model.weights, weight_permutation)

    sibling.location += parent.location
    parent.location = 0

    PruneIndexFromTree!(tree, prune_index)
end

# Graft tree while also adjusting the weight matrix appropriately
function graft_tree!(model::ModelState,
                     prune_index::Int,
                     graftpoint_index::Int,
                     parent_features::Array{Int64,1},
                     graftpoint_features::Array{Int64,1})

    tree = model.tree
    weight_indices = weight_index_pointers(tree)

    parent = tree.nodes[prune_index].parent
   
    graftpoint_start = weight_indices[graftpoint_index]
    parent_start = weight_indices[parent.index]

    graftnode = tree.nodes[graftpoint_index]
    graftpoint_end = graftpoint_start + graftnode.location - 1
    parent_end = parent_start + parent.location - 1

    weight_permutation = linspace(1, size(model.weights)[1], size(model.weights)[1])

    if graftpoint_start < parent_start
        for p = reverse(parent_features)
            insert(weight_permutation, parent_start, p)
        end

        del(weight_permutation, graftpoint_start:graftpoint_end)
        for p = reverse(graftpoint_features)
            insert(weight_permutation, graftpoint_start, p)
        end
    else
        del(weight_permutation, graftpoint_start:graftpoint_end)
        for p = reverse(graftpoint_features)
            insert(weight_permutation, graftpoint_start, p)
        end

        for p = reverse(parent_features)
            insert(weight_permutation, parent_start, p)
        end
    end

    permute_rows_and_cols(model.weights, weight_permutation)

    graftnode.location = length(graftpoint_features)
    parent.location = length(parent_features)

    InsertIndexIntoTree!(tree, prune_index, graftpoint_index)
end

