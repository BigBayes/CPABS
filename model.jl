load("tree.jl")

# Type for augmented weight matrix.  store a matrix, a set of pointers for 
# each node denoting features attributed to said node
type AugmentedMatrix
    matrix::Array{Float64,2}

    num_sets::Int
    feature_pointers::Array{Array{Int,1},1} 
    num_active_features::Array{Int,1}

    next_novel_index::Int
end

AugmentedMatrix(N::Int) = 
    AugmentedMatrix(zeros(10,10), N, [Int[] for x = 1:N], zeros(Int,N), 1)
function AugmentedMatrix(N::Int, W::Array{Float64,2}, feature_counts::Array{Int,1})

    feature_pointers = [Int[] for x = 1:N]
    next_index = 1 
    for i = 1:N
        append!(feature_pointers[i],linspace(next_index, next_index + feature_counts[i] - 1, feature_counts[i]))
        next_index += feature_counts[i]
    end
    AugmentedMatrix(copy(W), N, feature_pointers, feature_counts, size(W)[1]+1)
end

copy(aug_matrix::AugmentedMatrix) = AugmentedMatrix(copy(aug_matrix.matrix), 
                                                    aug_matrix.num_sets,
                                                    deepcopy(aug_matrix.feature_pointers),
                                                    copy(aug_matrix.num_active_features),
                                                    aug_matrix.next_novel_index)

ref(aug_matrix::AugmentedMatrix, rows, cols) = aug_matrix.matrix[rows,cols]

function assign(aug_matrix::AugmentedMatrix, array, rows, cols)
    aug_matrix.matrix[rows,cols] = array
end

type ModelState
    lambda::Float64
    gamma::Float64
    w_sigma::Float64
    b_sigma::Float64
    nu::Float64
    tree::Tree{Int64}
    weights::Array{Float64,2}
    augmented_weights::AugmentedMatrix
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

    # Inference Params
    rrj_jump_probabilities::Array{Float64} #assumes L \in {k-1,k,k+1}
    global_move_probability::Float64
    Z_sample_branch_prob::Float64

    debug::Bool
    verbose::Bool
end

copy(ms::ModelSpecification) = ModelSpecification(ms.use_pairwise, ms.use_parenthood,
                                   ms.use_childhood, ms.symmetric_W, ms.diagonal_W,
                                   copy(ms.rrj_jump_probabilities), 
                                   ms.global_move_probability,
                                   ms.Z_sample_branch_prob,
                                   ms.debug, ms.verbose)

# Data container
type DataState
    Ytrain::Array{Int64,2}
    Ytest::Array{Int64,2}
    X_r::Array{Float64,3}
    X_p::Array{Float64,2}
    X_c::Array{Float64,2}
end

copy(data::DataState) = DataState(copy(data.Ytrain), 
                                  copy(data.Ytest),
                                  copy(data.X_r),
                                  copy(data.X_p),
                                  copy(data.X_c)) 
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

    # adjust augmented matrix
    aug_W = model.augmented_weights
    num_features = aug_W.num_active_features[parent.index]
    move_features(aug_W,
                  linspace(1,num_features,num_features),
                  parent.index,
                  sibling.index) 

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

    # adjust augmented matrix
    aug_W = model.augmented_weights
    # augmented_weights feature offsets are relative, given ones are absolute
    aug_parent_features = parent_features - min([graftpoint_features, parent_features]) + 1
    move_features(aug_W, aug_parent_features, graftnode.index, parent.index)

    graftnode.state = length(graftpoint_features)
    parent.state = length(parent_features)

    InsertIndexIntoTree!(tree, prune_index, graftpoint_index)
end

function permute_rows_and_cols!(A, permutation)
    assert( length(permutation) == size(A)[1])
    A[:,:] = A[:,permutation]
    A[:,:] = A[permutation,:]
end

function delete_row_and_col(A, index)
    lenA = size(A)[1]
    valid_indices = [ x < index ? x : x+1 for x = 1:lenA-1]
    A[valid_indices,valid_indices]
end



function get_augmented_submatrix_indices(augmented_matrix::AugmentedMatrix,
                                         augmented_set_index::Int,
                                         num_augmented_features::Int)
    matrix_indices = Int[]
    feature_pointers = augmented_matrix.feature_pointers
    num_active_features = augmented_matrix.num_active_features

    for i = 1:length(feature_pointers)
        num_features = num_active_features[i]

        if i == augmented_set_index
            num_features += num_augmented_features
        end

        if num_features > length(feature_pointers[i])
            add_new_features(augmented_matrix,
                             i,
                             num_features - length(feature_pointers[i]) )
        end
        append!(matrix_indices, feature_pointers[i][1:num_features])
    end

    # always give augmented features last
#    i = augmented_set_index
#    if num_augmented_features  > 0
#        num_features = num_active_features[i] + num_augmented_features
#        if num_features > length(feature_pointers[i])
#            add_new_features(augmented_matrix,
#                             i,
#                             num_features - length(feature_pointers[i]) )
#        end
#        append!(matrix_indices, feature_pointers[i][num_active_features[i]+1:num_features])
#    end

    matrix_indices
end


function expand_matrix(augmented_matrix::AugmentedMatrix,
                       new_size::Int)
    old_size = size(augmented_matrix.matrix)[1]
    assert(new_size > old_size)
    A = zeros(new_size, new_size)
    A[1:old_size,1:old_size] = augmented_matrix.matrix
    augmented_matrix.matrix = A 
end

# adds new features (but not making them active)
function add_new_features(augmented_matrix::AugmentedMatrix,
                          augmented_set_index::Int,
                          num_features::Int)
    next_index = augmented_matrix.next_novel_index
    end_index = next_index + num_features - 1
    new_features = linspace(next_index,end_index,num_features)

    matrix_size = size(augmented_matrix.matrix)[1]
    if end_index > matrix_size 
        expand_matrix(augmented_matrix, end_index + 20)
    end

    append!(augmented_matrix.feature_pointers[augmented_set_index], new_features)
    augmented_matrix.next_novel_index = end_index + 1
end

function deactivate_feature(augmented_matrix::AugmentedMatrix,
                            augmented_set_index::Int,
                            feature_index::Int)
    feature_pointers = augmented_matrix.feature_pointers[augmented_set_index]
    num_features = augmented_matrix.num_active_features[augmented_set_index]
    el = feature_pointers[feature_index]

    del(feature_pointers, feature_index)
    insert(feature_pointers, num_features, el) #insert el as next augmented feature
    augmented_matrix.num_active_features[augmented_set_index] -= 1
end

function activate_feature(augmented_matrix::AugmentedMatrix,
                          augmented_set_index::Int)
    feature_pointers = augmented_matrix.feature_pointers[augmented_set_index]
    num_features = augmented_matrix.num_active_features
    num_features[augmented_set_index] += 1
    if num_features[augmented_set_index] > length(feature_pointers)
        add_new_features(augmented_matrix, augmented_set_index, 1)
    end
end

# for moving features when sampling psi changes feature assignments
function move_features(augmented_matrix,
                       feature_indices::Array{Int,1},
                       source_index::Int,
                       target_index::Int)
    source_pointers = augmented_matrix.feature_pointers[source_index]
    target_pointers = augmented_matrix.feature_pointers[target_index]

    num_features = augmented_matrix.num_active_features
    els = Int[]

    rev_indices = reverse(sort(feature_indices))
    for feature_index = rev_indices
        el = source_pointers[feature_index]
        push(els, el)
        del(source_pointers, feature_index)
        num_features[source_index] -= 1
    end

    for el = reverse(els)
        numf = num_features[target_index]
        insert(target_pointers, numf+1, el)
        num_features[target_index] += 1
    end
end

# Utilities for saving/restoring models

function model2array(model::ModelState)
    (Z, U, inds) = tree2array(model.tree, model.gamma)
    _2Nm1 = length(inds)
    N::Int = (_2Nm1+1)/2
    feature_indices = [zeros(Int, 0) for i = 1:2N-1]
    W_index_pointers = weight_index_pointers(model.tree)

    for i = 1:2N-1
        start_index = W_index_pointers[i]
        u = model.tree.nodes[i].state
        end_index = start_index + u - 1

        feature_indices[i] = linspace(start_index,end_index,u) 
    end

    permuted_features = feature_indices[inds]
   
    permutation = Int64[]


    for i = 1:2N-1
        append!(permutation, permuted_features[i])
    end

    W = copy(model.weights)
    permute_rows_and_cols!(W, permutation)

    (Z, U, W)
end
