load("model.jl")
load("tree.jl")
load("probabiltiy_util.jl")

#Global pdf
function full_pdf(model::ModelState,
                  Y::Array{Int64,2},
                  X_r::Array{Float64,3},
                  X_p::Array{Float64,2},
                  X_c::Array{Float64,2})
    (N,N) = size(Y)
    prior(model)+likelihood(model,Y,X_r,X_p,X_c,1:N)
end

function prior(model::ModelState)
    total_prob = sum(normal_logpdf(model.weights))
    total_prob += sum(normal_logpdf(model.beta))
    total_prob += sum(normal_logpdf(model.beta_p))
    total_prob += sum(normal_logpdf(model.beta_r))
    total_prob += sum(normal_logpdf(model.a))
    total_prob += sum(normal_logpdf(model.b))
    total_prob += normal_logpdf(model.c)

    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N = (_2Nm1+1)/2
    total_prob += sum(log(1:N))
    for i = 1:length(tree.nodes)
        total_prob -= log(tree.nodes[i].num_leaves)
        segment_length = (1-gamma)*gamma^(tree.nodes[i].num_ancestors-1)
        total_prob += poisson_logpdf(tree.nodes[i].location, lambda * segment_length)
    end
    total_prob
end

function likelihood(model::ModelState,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2},
                    indices::Array{Int64,1})
    (N,N) = size(Y)
    Z = ConstructZ(tree)
    W = model.weights
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_r
 
    total_prob = 0 
    for i = indices
        for j = 1:N
            logit_arg = Z[i,:] * W * Z[j,:]' + 
                        beta' * X_r[i,j] + beta_p' * X_p[i] + beta_c' * X_c[j] + 
                        model.a[i] + model.b[j] + model.c
            if Y[i,j] == 1
                total_prob += -log(1+exp(-logit_arg))
            elseif Y[i,j] == 0
                total_prob += -logit_arg - log(1+exp(-logit_arg))
            end
        end
    end
    total_prob
end

#Local pdfs for sampling

#pdf for updating psi, the tree structure, assumes tree is already pruned at prune_index
function infsites_logpdf(model::ModelState,
                         pruned_index::Int64)
    (N,N) = size(Y)
    tree = model.tree

    #psi_probs = prior_path(tree, pruned_index, path_leaf_index)
    psi_probs = prior_tree(tree, pruned_index)

    descendant_mutation_probs = -Inf * ones(2N - 1)
    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    total_probs = Float64[]
    tree_states = Array(Tuple,0)

    i = 1
    while contains(subtree_indices, i)
        i += 1
    end
    root = FindRoot(tree, i)

    indices = GetLeafToRootOrdering(tree, root)

    subtree_probs = -Inf * ones(2N - 1)
    for i = indices
        graft_node = tree.nodes[i]
        for j = subtree_indices
            cur = tree.nodes[j]
            if j <= N
                poisson_mean = model.lambda * gamma^(cur.num_ancestors + graft_node.num_ancestors)
            else
                poisson_mean = model.lambda * (1 - gamma) * gamma^(cur.num_ancestors + graft_node.num_ancestors - 1)
            end
            subtree_probs[i] += poisson_logpdf(cur.location, poisson_mean)
        end
    end

    for i = indices
        cur = tree.nodes[i]

        if i == 1 #leaf node
            poisson_mean_before = model.lambda * gamma^(cur.num_ancestors)
        else
            poisson_mean_before = model.lambda * (1 - gamma) * gamma^(cur.num_ancestors - 1)
        end

        poisson_mean_after = poisson_mean_before * gamma
        descendant_mutation_probs[i] += poisson_logpdf(cur.location, poisson_mean_after) -
                                        poisson_logpdf(cur.location, poisson_mean_before)

        child_mutation_contributions = 0.0
        if i > 1
            for j = 1:2
                child = cur.children[j]
                if child != Nil()
                    assert( descendant_mutation_probs[child.index] != -Inf)
                    child_mutation_contributions += descendant_mutation_probs[child.index]
                end
            end
        end

        descendant_mutation_probs[i] += child_mutation_contributions

        (U, V) = all_splits(linspace(1, cur.location, cur.location))

        for j = 1:length(U)
            u = U[j]
            v = V[j]

            prob = poisson_logpdf(length(u), poisson_mean_after) +
                   poisson_logpdf(length(v), poisson_mean_before) -
                   poisson_logpdf(cur.location, poisson_mean_before)

            push(total_probs, prob + child_mutation_contributions + psi_probs[i] + subtree_probs[i])

            # (node above which to attach, features moved below graft point, features moved above)
            push(tree_states, (i, u, v) )                

        end
    end

    (total_probs, tree_states)
end

# Likelihood for mcmc moves for psi
# Assumes features are moved around appropriately after pruning the tree (ie the
# features above the pruned parent are moved to the unpruned child
function psi_likelihood_logpdf(model::ModelState,
                               pruned_index::Int64,
                               Y::Array{Int64,2},
                               X_r::Array{Float64,3},
                               X_p::Array{Float64,2},
                               X_c::Array{Float64,2})
    (N,N) = size(Y)
    tree = model.tree

    weight_indices = weight_index_pointers(tree)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)
    subtree_leaves = subtree_indices[find(subtree_indices <= N)]

    i = 1
    while contains(subtree_indices, i)
        i += 1
    end
    root = FindRoot(tree, i)

    indices = GetLeafToRootOrdering(tree, root.index)
    leaves = GetLeaves(tree, root.index)

    tree_states = Array(Tuple,0)

    observed_parenthood_effects = zeros(Float64, (length(subtree_leaves), N))
    observed_childhood_effects = zeros(Float64, (length(subtree_leaves), N))

    subtree_features = fill(Array(Int64,0), length(subtree_leaves))

    for i = subtree_leaves
        for j = 1:N
            observed_parenthood_effects[i,j] = 
                model.beta' * squeeze(X_r[i,j,:]) + 
                model.beta_p' * squeeze(X_p[i,:]) + 
                model.beta_c' * squeeze(X_c[j,:]) +
                model.a[i] + model.b[j] + model.c

            observed_childhood_effects[i,j] = 
                model.beta' * squeeze(X_r[j,i,:]) + 
                model.beta_p' * squeeze(X_p[j,:]) + 
                model.beta_c' * squeeze(X_c[i,:]) +
                model.a[j] + model.b[i] + model.c
        end
    end

    for i = subtree_indices
        cur = tree.nodes[i]
        features_start = weight_indices[i]
        features_end = features_start + cur.location - 1

        cur_features = linspace(features_start, features_end, cur.location)

        # cycle through all subtree leaves, find if the current node is its ancestor
        # and if so add the current node's features to the leaf's features for 
        # likelihood computations.  Can certainly be optimized
        for l = 1:length(subtree_leaves)
            ancestor = tree.nodes[subtree_leaves[l]]
            while ancestor != Nil()
                if ancestor.index == i
                    append!(subtree_features[l], cur_features)
                end
                ancestor = ancestor.parent
            end
        end
    end 

    latent_parenthood_effects = zeros(Float64, (length(subtree_leaves, N)))
    latent_childhood_effects = zeros(Float64, (length(subtree_leaves, N)))


    Z = constructZ(tree)

    leaf_features = fill(Array(Int64,0), length(1:N))
    for l = leaves
        leaf_features[l] = find(Z[l,:] > 0)
    end


    likelihoods = Float64[]
    for i = indices
        cur = tree.nodes[i]

        ancestors = GetPath(tree, i)
        features = Int64[]

        for k = 1:length(ancestors)
            acur = ancestors[k]
            features_start = weight_indices[acur.index]
            features_end = features_start + acur.location - 1

            append!(features, linspace(features_start, features_end, acur.location))
        end
      
        features_start =  weight_indices[i]
        features_end = weight_indices[i] + cur.location - 1

        (V, T) = all_splits(linspace(features_start, features_end, cur.location))
        for m = 1:length(T)
            # t is the set of features to be included above the graft point,
            # so they are included in the likelihood computation for the 
            # leaves of the pruned subtree
            t = T[m]

            likelihood = 0.0

            # TODO Optimize this if necessary (eg the appends in the second inner for loop can be removed)
            for j1 = 1:length(subtree_leaves)
                l1 = subtree_leaves[j1]
                total_features = Int64[]
                append!(total_features, features)
                append!(total_features, t)
                append!(total_features, subtree_features[j1])
                for l2 = leaves
                    latent_effect = sum(W[total_features, leaf_features[l2]])
                    likelihood += log_logit(latent_effect + observed_parenthood_effects[l1,l2] + observed_childhood_effects[l1,l2], Y[l1,l2])
                end

                for j2 = 1:length(subtree_leaves)
                    l2 = subtree_leaves[j2]
                    subtree_leaf_features = Int64[]
                    append!(subtree_leaf_features, features)
                    append!(subtree_leaf_features, t)
                    append!(subtree_leaf_features, subtree_features[j2])
                    latent_effect = sum(W[total_features, subtree_leaf_features])
                    likelihood += log_logit(latent_effect + observed_parenthood_effects[l1,l2] + observed_childhood_effects[l1,l2], Y[l1,l2])
                end
            end
            push(likelihoods, likelihood)
            push(tree_states, (i, V[m], t))
        end       

    end
    (likelihoods, tree_states)
end

# log \prod M_i for all valid grafting points
function prior_tree(tree::Tree,
                    pruned_index::Int)
    N = (length(tree.nodes) + 1) / 2

    # -log \prod M_i for all points below current one, both original and 
    # for grafted trees where the subtree was grafted below the current node
    original_probs = -Inf*ones(2N - 1)
    new_probs = -Inf*ones(2N - 1)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    i = 1
    while contains(subtree_indices, i)
        i += 1
    end
    root = FindRoot(tree, i)

    # indices will contain a leaves to root ordering of nodes in the pruned tree
    indices = GetLeafToRootOrdering(tree, root.index)

    for i = indices
        cur = tree.nodes[i]
        original_probs[i] = -log(max(cur.num_leaves - 1, 1))
        new_probs[i] = -log(cur.num_leaves + tree.nodes[pruned_index].num_leaves - 1)
        if i > N
            for j = 1:2
                child_original_prob = original_probs[cur.children[j].index]
                child_new_prob = new_probs[cur.children[j].index]
                if child_original_prob == -Inf || child_new_prob == -Inf
                    assert(contains(subtree_indices, cur.children[j].index))
                end
                original_probs[i] += child_original_prob
                new_probs[i] += child_new_prob
            end
        end
    end

    new_probs - original_probs
end

# log \prod M_i for a path starting at leaf node_index
function prior_path(tree::Tree,
                    pruned_ind::Int,
                    node_index::Int)
    num_leaves = GetPathNumDescendants(tree, node_index)
    proposed_num_leaves = ProposedNumDescendants(tree,pruned_ind,node_index)
    original_partial_values = cumsum(log(num_leaves-1))
    proposed_partial_values = cumsum(log(proposed_num_leaves-1))
    -(original_partial_values[end]-original_partial_values+proposed_partial_values)
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

# Utility Functions

# constructs all possible splits of v into 2 sets
function all_splits(v::Array{Int64})
    # Construct power set U, and for each u \in U, construct w = v \ u
    # Then the set of (u,w) will be all possible splits of v

    assert(length(v) < 20)

    U = Vector{Int64}[]
    
    # Construct power set U
    for i = 0:2^(length(v))-1
        x = i
        u = Int64[]
        for j = 1:length(v)
            if bool(x & 1)
                push(u, v[j])
            end
            x = x >> 1
        end
        push(U, u)
    end
    (U,reverse(U))
end

# Pointers to the beginning of a node's weights in model.weights
function weight_index_pointers(tree::Tree)
    n = length(tree.nodes) 

    U = zeros(Int64, n)
    for i = 1:n
        U[i] = tree.nodes[i].location
    end

    # Index pointing to starting location in W for the features found
    # above particular nodes
    weight_indices = cumsum(U)
    weight_indices[2:end] = weight_indices[1:end-1]
    weight_indices[1] = 0
    weight_indices += 1

    weight_indices
end

function permute_rows_and_cols!(A, permutation)
    A[:,:] = A[:,permutation]
    A[:,:] = A[permutation,:]
end

