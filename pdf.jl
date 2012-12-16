load("model.jl")
load("tree.jl")
load("probability_util.jl")
load("profile.jl")
#Global pdf
function full_pdf(model::ModelState,
                  model_spec::ModelSpecification,
                  Y::Array{Int64,2},
                  X_r::Array{Float64,3},
                  X_p::Array{Float64,2},
                  X_c::Array{Float64,2})
    (N,N) = size(Y)
    prior(model)+likelihood(model,model_spec,Y,X_r,X_p,X_c,linspace(1,N,N))
end

function prior(model::ModelState)
    w_sigma = model.w_sigma
    b_sigma = 1.0
    total_prob = sum(normal_logpdf(model.weights, w_sigma ))
    total_prob += sum(normal_logpdf(model.beta, b_sigma))
    total_prob += sum(normal_logpdf(model.beta_p, b_sigma))
    total_prob += sum(normal_logpdf(model.beta_c, b_sigma))
    total_prob += sum(normal_logpdf(model.a, b_sigma))
    total_prob += sum(normal_logpdf(model.b, b_sigma))
    total_prob += normal_logpdf(model.c, b_sigma)

    gam = model.gamma
    lambda = model.lambda
    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2
    total_prob += sum(log(1:N))
    for i = 1:length(tree.nodes)
        total_prob -= log(tree.nodes[i].num_leaves)
        segment_length = (1-gam)*gam^(tree.nodes[i].num_ancestors-1)
        total_prob += poisson_logpdf(tree.nodes[i].state, lambda * segment_length)
    end
    total_prob
end

function likelihood(model::ModelState,
                    model_spec::ModelSpecification,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2},
                    indices::Array{Int64,1})
    (N,N) = size(Y)
    tree = model.tree
    Z = ConstructZ(tree)
    W = model.weights
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c
 
    total_prob = 0.0
    for i = indices
        for j = 1:N

#            if i == 1 && mod(j, 20) == 0
#                println("likelihood latent effect: ", squeeze(Z[i,:] * W * Z[j,:]')[1] )
#            end
            oe = compute_observed_effects(model, model_spec, i, j, X_r, X_p, X_c)
            logit_arg = (Z[i,:] * W * Z[j,:]' + oe)[1]
                        

            total_prob += log_logit(logit_arg, Y[i,j])
        end
    end
    total_prob
end

# Local pdfs for sampling

###################################
###### pdfs for psi updates #######
###################################

#@profile begin
#pdf for updating psi, the tree structure, assumes tree is already pruned at prune_index
function psi_infsites_logpdf(model::ModelState,
                             pruned_index::Int64,
                             path::Array{Int64,1})
    (N,N) = size(Y)
    tree = model.tree
    gam = model.gamma
    lambda = model.lambda
    weight_indices = weight_index_pointers(tree)

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

    indices = GetLeafToRootOrdering(tree, root.index)

    subtree_probs = -Inf * ones(2N - 1)
    for i = indices
        graft_node = tree.nodes[i]
        if length(subtree_indices) > 0
            subtree_probs[i] = 0.0
        end
        for j = subtree_indices
            cur = tree.nodes[j]
            if j <= N
                poisson_mean = model.lambda * gam^(cur.num_ancestors + graft_node.num_ancestors)
            else
                poisson_mean = model.lambda * (1 - gam) * gam^(cur.num_ancestors + graft_node.num_ancestors - 1)
            end
            subtree_probs[i] += poisson_logpdf(cur.state, poisson_mean)
        end
    end

    for i = indices
        cur = tree.nodes[i]

        if i == 1 #leaf node
            poisson_mean_before = model.lambda * gam^(cur.num_ancestors)
        else
            poisson_mean_before = model.lambda * (1 - gam) * gam^(cur.num_ancestors - 1)
        end

        poisson_mean_after = poisson_mean_before * gam
        descendant_mutation_probs[i] = poisson_logpdf(cur.state, poisson_mean_after) -
                                       poisson_logpdf(cur.state, poisson_mean_before)

        child_mutation_contributions = 0.0
        if i > 1
            for j = 1:2
                child = cur.children[j]
                if child != Nil()
                    if descendant_mutation_probs[child.index] == -Inf
                        println("i: ", i)
                        println("indices: ", indices)
                        cur_index = find(indices .== cur.index)[1]
                        child_index = find(indices .== child.index)[1]
                        println("(cur, child): ", (cur.index, child.index)) 
                        assert(false)
                    end
                    child_mutation_contributions += descendant_mutation_probs[child.index]
                end
            end
        end

        descendant_mutation_probs[i] += child_mutation_contributions

        features_start = weight_indices[i]
        features_end = weight_indices[i] + cur.state - 1

        if contains(path, i)
            (U, V) = all_splits(linspace(features_start, features_end, cur.state))
            for j = 1:length(U)
                u = U[j]
                v = V[j]

                prob = poisson_logpdf(length(u), poisson_mean_after) +
                       poisson_logpdf(length(v), poisson_mean_before) -
                       poisson_logpdf(cur.state, poisson_mean_before)

                push(total_probs, prob + child_mutation_contributions + psi_probs[i] + subtree_probs[i])
    #            if true
    #                println("prob ", prob)
    #                println("child contributions ", child_mutation_contributions)
    #                println("psi prob ", psi_probs[i])
    #                println("subtree_probs ", subtree_probs[i])
    #            end
                # (node above which to attach, features moved below graft point, features moved above)
                push(tree_states, (i, u, v) )                
            end
        end
    end

    (total_probs, tree_states)
end

# Likelihood for mcmc moves for psi
# Assumes features are moved around appropriately after pruning the tree (ie the
# features above the pruned parent are moved to the unpruned child
function psi_likelihood_logpdf(model::ModelState,
                               model_spec::ModelSpecification,
                               pruned_index::Int64,
                               path::Array{Int64,1},
                               Y::Array{Int64,2},
                               X_r::Array{Float64,3},
                               X_p::Array{Float64,2},
                               X_c::Array{Float64,2})
    (N,N) = size(Y)
    tree = model.tree
    W = model.weights

    weight_indices = weight_index_pointers(tree)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)
    subtree_indices_array = [k for k in subtree_indices]
    subtree_leaves = subtree_indices_array[find(subtree_indices_array .<= N)]

    #println("subtree_leaves: ", subtree_leaves)

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

    #subtree_features = fill(Array(Int64,0), length(subtree_leaves))
    subtree_features = [Array(Int64,0) for x = 1:length(subtree_leaves)]


    for i = 1:length(subtree_leaves)
        l1 = subtree_leaves[i]
        for l2 = 1:N
            observed_parenthood_effects[i,l2] = compute_observed_effects(model, model_spec, l1, l2, X_r, X_p, X_c)
            observed_childhood_effects[i,l2] = compute_observed_effects(model, model_spec, l2, l1, X_r, X_p, X_c)
        end
    end

    for i = subtree_indices
        cur = tree.nodes[i]
        features_start = weight_indices[i]
        features_end = features_start + cur.state - 1

        cur_features = linspace(features_start, features_end, cur.state)

        # cycle through all subtree leaves, find if the current node is its ancestor
        # and if so add the current node's features to the leaf's features for 
        # likelihood computations.  Could be optimized
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


    Z = ConstructZ(tree)

    leaf_features = [Array(Int64,0) for x = 1:N] #fill(Array(Int64,0), N) this is an array of references!!! 

    nzI, nzJ = findn(Z .> 0)

#    for l = leaves
#        leaf_features[l] = find(Z[l,:] .> 0)
#    end

    for i = 1:length(nzI)
        push(leaf_features[nzI[i]], nzJ[i])
    end

    latent_effects = zeros(Float64, (length(path), length(subtree_leaves), N))
    #target_subtree_features = fill(Array(Int64,0), (length(path), length(subtree_leaves)))
    target_subtree_features = [Array(Int64,0) for x = 1:length(path), y = length(subtree_leaves)]


    likelihoods = Float64[]
    rpath = reverse(path)
    for path_index = 1:length(path)
        i = rpath[path_index]
        cur = tree.nodes[i]

        if path_index > 1
            parent = cur.parent
            parent_path_index = find(rpath .== parent.index)[1]
        end

        #latent_effects = zeros(Float64, (length(subtree_leaves), N))
        #target_subtree_features = fill(Array(Int64,0), length(subtree_leaves))
        for j1 = 1:length(subtree_leaves)

            if path_index == 1

                assert(cur.parent == Nil())
                append!(target_subtree_features[path_index,j1], subtree_features[j1])

                for j2 = 1:length(leaves)
                    l2 = leaves[j2]
                    assert(subtree_leaves[j1] != l2)
                    latent_effects[path_index,j1,l2] = sum(W[target_subtree_features[path_index,j1], leaf_features[l2]])
                end
            else
                features_start = weight_indices[parent.index]
                features_end = features_start + parent.state - 1

                new_features = linspace(features_start, features_end, parent.state)
                target_subtree_features[path_index,j1] = copy(target_subtree_features[parent_path_index,j1])
                append!(target_subtree_features[path_index,j1], new_features)

                for j2 = 1:length(leaves)
                    l2 = leaves[j2]
                    latent_effects[path_index,j1,l2] = sum(W[new_features, leaf_features[l2]]) +
                                                       latent_effects[parent_path_index,j1,l2]
                end
            end

        end

        for j1 = 1:length(subtree_leaves)
            for j2 = 1:length(subtree_leaves)
                l2 = subtree_leaves[j2]
                latent_effects[path_index,j1,l2] = sum(W[target_subtree_features[path_index,j1], target_subtree_features[path_index,j2]])
            end
        end

        features_start = weight_indices[i]
        features_end = weight_indices[i] + cur.state - 1

        (V, T) = all_splits(linspace(features_start, features_end, cur.state))

        #do everything in reverse
        rT = reverse(T)
        rV = reverse(V)
        for m = 1:length(T)
            # t is the set of features to be included above the graft point,
            # so they are included in the likelihood computation for the 
            # leaves of the pruned subtree
            t = rT[m]

            likelihood = 0.0

            # TODO Optimize this if necessary (eg the appends in the second inner for loop could be removed)
            for j1 = 1:length(subtree_leaves)
                l1 = subtree_leaves[j1]
                #total_features = Int64[]
                #append!(total_features, features)
                #append!(total_features, t)
                #append!(total_features, subtree_features[j1])
                for l2 = leaves
                    assert(l1 != l2)
                    latent_effect = sum(W[t, leaf_features[l2]]) + latent_effects[path_index,j1,l2]

#                    if l1 == 1 && mod(l2,20) == 0
#                        println("latent effect: ", latent_effect)
#                        println("leaf features: ", leaf_features[l2])
#                    end

                    likelihood += log_logit(latent_effect + observed_parenthood_effects[j1,l2], Y[l1,l2])
                    likelihood += log_logit(latent_effect + observed_childhood_effects[j1,l2], Y[l2,l1])
                    if likelihood == -Inf
                        println("leaves")
                        println(latent_effect,",",observed_parenthood_effects[j1,l2],",",observed_childhood_effects[j1,l2])
                    end
                end

                for j2 = 1:length(subtree_leaves)
                    l2 = subtree_leaves[j2]
                    subtree_leaf_features1 = Int64[]
                    subtree_leaf_features2 = Int64[]
                    if l2 != l1
                        append!(subtree_leaf_features1, target_subtree_features[path_index,j1])
                        append!(subtree_leaf_features2, target_subtree_features[path_index,j2])
                    end
                    latent_effect = sum(W[t, subtree_leaf_features2]) + sum(W[t, subtree_leaf_features1]) + sum(W[t,t]) + latent_effects[path_index,j1,l2]
                    likelihood += log_logit(latent_effect + observed_parenthood_effects[j1,l2], Y[l1,l2])
                    #likelihood += log_logit(latent_effect + observed_childhood_effects[j1,l2], Y[l2,l1])
                    if likelihood == -Inf
                        println("subtree leaves")
                        println(latent_effect,",",observed_parenthood_effects[j1,l2],",",observed_childhood_effects[j1,l2])
                    end
                
                end
    
            end
            insert(likelihoods, 1, likelihood)
            insert(tree_states, 1, (i, rV[m], t))
        end       

    end

    (likelihoods, tree_states)
end

# log \prod M_i for all valid grafting points
function prior_tree(tree::Tree,
                    pruned_index::Int)
    N::Int = (length(tree.nodes) + 1) / 2

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

###################################
###### pdfs for W updates #########
###################################

# pdf for updating one element of W
function W_local_logpdf(model::ModelState,
                        Y::Array{Int64, 2},
                        relevant_pairs::Array{Int64, 2},
                        latent_effects::Array{Float64, 2},
                        observed_effects::Array{Float64, 2},
                        w_old::Float64,
                        w_new::Float64)
    sigma = model.w_sigma 

    (_, npairs) = size(relevant_pairs)
  
    logprob = 0.0 
    for p = 1:npairs
        i = relevant_pairs[1,p]
        j = relevant_pairs[2,p]

        le = latent_effects[i,j]
        oe = observed_effects[i,j]

        le += w_new - w_old

        logprob += log_logit(le + oe, Y[i,j])
    end 

    return logprob + normal_logpdf(w_new, sigma)
end

###################################
###### pdfs for vardim updates ####
###################################

function vardim_local_logpdf(model::ModelState,
                             Y::Array{Int64, 2},
                             relevant_pairs::Array{Int64, 2},
                             latent_effects::Array{Array{Float64, 2}, 1},
                             observed_effects::Array{Float64, 2},
                             mixture_component_index::Int64,
                             u::Int64,
                             node_index::Int64,
                             effective_lambda::Float64,
                             w_old::Float64,
                             w_new::Float64,
                             w_is_auxiliary::Bool)

    (K,K) = size(model.weights)
    sigma = model.w_sigma
    if u == model.tree.nodes[node_index].state - 1
        if w_is_auxiliary
            logprob = W_local_logpdf(model,Y,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            # W_local_logpdf includes the prior term for w_old as if it were normal,
            # need to subtract it out here
            logprob += t_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,Y,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    elseif u == model.tree.nodes[node_index].state
        if w_is_auxiliary
            logprob = W_local_logpdf(model,Y,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            logprob += t_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,Y,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    else
        logprob = W_local_logpdf(model,Y,relevant_pairs, latent_effects[mixture_component_index],
                  observed_effects, w_old, w_new)
    end

    logprob + poisson_logpdf(u, effective_lambda)
end

function vardim_local_splits(model::ModelState,
                             Y::Array{Int64, 2},
                             relevant_pairs::Array{Int64, 2},
                             latent_effects::Array{Array{Float64, 2}, 1},
                             observed_effects::Array{Float64, 2},
                             u::Array{Int64, 1},
                             weight_index_pointers::Array{Int64, 1},
                             node_index::Int64,
                             effective_lambda::Float64,
                             constant_terms::Array{Float64,1},
                             w_old::Float64,
                             w_new::Float64,
                             w_is_auxiliary::Array{Bool, 1})

    num_components = length(latent_effects)
    logprobs = zeros(num_components)

    for mixture_component = 1:num_components
        logprob = vardim_local_logpdf(model, Y, relevant_pairs, latent_effects,
                                      observed_effects, mixture_component, 
                                      u[mixture_component], node_index,
                                      effective_lambda, w_old, w_new,
                                      w_is_auxiliary[mixture_component]) 
        logprobs[mixture_component] = logprob
    end
    logprobs + constant_terms
end

function vardim_local_sum(model::ModelState,
                          Y::Array{Int64, 2},
                          relevant_pairs::Array{Int64, 2},
                          latent_effects::Array{Array{Float64, 2}, 1},
                          observed_effects::Array{Float64, 2},
                          u::Array{Int64, 1},
                          weight_index_pointers::Array{Int64, 1},
                          node_index::Int64,
                          effective_lambda::Float64,
                          constant_terms::Array{Float64,1},
                          w_old::Float64,
                          w_new::Float64,
                          w_is_auxiliary::Array{Bool, 1})
    logprobs = vardim_local_splits(model, Y, relevant_pairs, latent_effects,
                                   observed_effects, u, weight_index_pointers,
                                   node_index, effective_lambda, constant_terms, 
                                   w_old, w_new, w_is_auxiliary)
    logsumexp(logprobs)
end

function vardim_logpdf(model::ModelState,
                       Y::Array{Int64, 2},
                       latent_effects::Array{Float64, 2},
                       observed_effects::Array{Float64, 2},
                       u::Int64,
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64,
                       w_is_auxiliary::Array{Bool, 2})
    (N,N) = size(Y)

    logprob = 0.0
    for i = 1:N
        for j = 1:N
            le = latent_effects[i,j]
            oe = observed_effects[i,j]
            logprob += log_logit(le + oe, Y[i,j])   
        end
    end


    W = model.weights
    nu = model.nu
    for k = find(w_is_auxiliary)
        logprob += t_logpdf(W[k], nu)
    end

    sigma = model.w_sigma
    for k = find(!w_is_auxiliary)
        logprob += normal_logpdf(W[k], sigma)
    end

    logprob + poisson_logpdf(u, effective_lambda)
end
function vardim_splits(model::ModelState,
                       Y::Array{Int64, 2},
                       latent_effects::Array{Array{Float64, 2}, 1},
                       observed_effects::Array{Float64, 2},
                       U::Array{Int64, 1},
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64,
                       w_is_auxiliary::Array{Array{Bool, 2}, 1})
    logprobs = zeros(length(U))
    for u_ind = 1:length(U)
        u = U[u_ind]
        logprobs[u_ind] = vardim_logpdf(model,Y, latent_effects[u_ind], observed_effects,
                                        u, weight_index_pointers, node_index,
                                        effective_lambda, w_is_auxiliary[u_ind])
    end
    logprobs
end
# Utility Functions

# Find all pairs i,j such that Z[i,k1]*Z[j,k2] = 1
function compute_relevant_pairs(Z, #sparse or full binary array
                                k1::Int64,
                                k2::Int64)
    ZZ = Z[:,k1]*Z[:,k2]'
    #ZZ = ZZ + ZZ' - ZZ .* ZZ'
    (I, J) = findn(ZZ)
    [I', J']
end

#function compute_relevant_pairs(Z,
#                                k1::Int64,
#                                k2::Int64)
#    k1_I = find(Z[:,k1])
#    k2_I = find(Z[:,k2])
#
#    I = Int64[]
#    J = Int64[]
#
#    for i = k1_I
#        for j = k2_I
#            push(I,i)
#            push(J,j)
#        end
#    end
#
#    [I', J']
#
#end

# Doesn't quite handle the symmetric case yet (eg the a[i] and b[j])
function compute_observed_effects(model::ModelState,
                                  model_spec::ModelSpecification,
                                  i::Int64,
                                  j::Int64,
                                  X_r::Array{Float64,3},
                                  X_p::Array{Float64,2},
                                  X_c::Array{Float64,2})

    observed_effect = model.c

    if model_spec.use_pairwise
        observed_effect += model.beta' * X_r[i,j,:]
    end

    if model_spec.use_parenthood
        observed_effect += model.beta_p' * X_p[i,:] + model.a[i]
    end

    if model_spec.use_childhood
        observed_effect += model.beta_c' * X_c[j,:] + model.b[j]
    end                                  

    observed_effect
end

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
        U[i] = tree.nodes[i].state
    end

    # Index pointing to starting location in W for the features found
    # above particular nodes
    weight_indices = cumsum(U)
    weight_indices[2:end] = weight_indices[1:end-1]
    weight_indices[1] = 0
    weight_indices += 1

    weight_indices
end

#end #profile

