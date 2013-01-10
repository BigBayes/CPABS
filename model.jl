load("tree.jl")

type ModelState
    lambda::Float64
    gamma::Float64
    w_sigma::Float64
    b_sigma::Float64
    nu::Float64
    tree::Tree{Int64}
    weights::Array{Float64,2}
    augmented_weights::Array{Float64,2}
    beta::Array{Float64,1}
    beta_p::Array{Float64,1}
    beta_c::Array{Float64,1}
    a::Array{Float64,1}
    b::Array{Float64,1}
    c::Float64
end

function copy(model::ModelState)
    ModelState(model.lambda, model.gamma, model.w_sigma, model.b_sigma, 
               model.nu, copy(model.tree), copy(model.weights),
               copy(model.augmented_weights), copy(model.beta), 
               copy(model.beta_p), copy(model.beta_c), 
               copy(model.a), copy(model.b), model.c)
end

# Model and inference specifications
type ModelSpecification
    use_pairwise::Bool
    use_parenthood::Bool
    use_childhood::Bool
    symmetric_W::Bool
    diagonal_W::Bool

    global_move_probability::Float64
    Z_sample_branch_prob::Float64

    debug::Bool
    verbose::Bool
end

# Data container
type DataState
    Ytrain::Array{Int64,2}
    Ytest::Array{Int64,2}
    X_r::Array{Float64,3}
    X_p::Array{Float64,2}
    X_c::Array{Float64,2}
end

# Tree/Model Utility Functions

# Prune tree while also adjusting the weight matrix appropriately
function prune_tree!(model::ModelState,
                     prune_index::Int)
    tree = model.tree
    weight_indices = weight_index_pointers(tree)

    parent = tree.nodes[prune_index].parent

    if parent.children[1].index == prune_index
        sibling = parent.children[2]
    else
        sibling = parent.children[1]
    end


    sibling_start = weight_indices[sibling.index]
    parent_start = weight_indices[parent.index]

    sibling_end = sibling_start + sibling.state - 1
    parent_end = parent_start + parent.state - 1

    weight_permutation = linspace(1,size(model.weights)[1],size(model.weights)[1])


    if sibling_start < parent_start
        del(weight_permutation, parent_start:parent_end)
    
        for p = reverse(parent_start:parent_end)
            insert(weight_permutation, sibling_end + 1, p)
        end
    else
        for p = reverse(parent_start:parent_end)
            insert(weight_permutation, sibling_end + 1, p)
        end

        del(weight_permutation, parent_start:parent_end)
    end

    permute_rows_and_cols!(model.weights, weight_permutation)

    sibling.state += parent.state
    parent.state = 0

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
    graftpoint_end = graftpoint_start + graftnode.state - 1
    parent_end = parent_start + parent.state - 1

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

    permute_rows_and_cols!(model.weights, weight_permutation)

    graftnode.state = length(graftpoint_features)
    parent.state = length(parent_features)

    InsertIndexIntoTree!(tree, prune_index, graftpoint_index)
end

function permute_rows_and_cols!(A, permutation)
    assert( length(permutation) == size(A)[1])
    A[:,:] = A[:,permutation]
    A[:,:] = A[permutation,:]
end

