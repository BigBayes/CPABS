require("model.jl")
require("tree.jl")
require("probability_util.jl")

using NumericExtensions
import NumericExtensions.evaluate, NumericExtensions.result_type

#require("profile.jl")
#Global pdf
function full_pdf(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState)
    YY = data.Ytrain
    (N,N) = size(YY[1])
    prior(model,model_spec)+likelihood(model,model_spec,data,[1:N])[1]
end

function prior(model::ModelState,
               model_spec::ModelSpecification)
    w_sigma = model.w_sigma
    b_sigma = model.b_sigma

    
    if model_spec.positive_W
        W = exp(model.weights)
        W_pdf = exp_logpdf
        W_jacobian = (W) -> W
    else
        W = model.weights
        W_pdf = normal_logpdf
        W_jacobian = (W) -> 0.0
    end
 
    if model_spec.diagonal_W && length(model.weights) > 0
        total_prob = sum(W_pdf(diag(W), w_sigma ))
    elseif model_spec.symmetric_W
        total_prob = sum(tril(W_pdf(model.weights, w_sigma )))
    else
        total_prob = sum(W_pdf(model.weights, w_sigma ))
    end

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

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    inf_term = 0.0
    psi_term = 0.0

    mu = ones(2N - 1)
    t = ones(2N - 1)
    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            mu[i] = mu[parent.index]*cur_mu_prop
            t[i] = t[parent.index]*(parent.rhot*cur_mu_prop)^gam
        end
    end

    t[1:N] = 0.0

    for i = 1:length(tree.nodes)
        cur = tree.nodes[i]
        parent = cur.parent
        if i > N
            left_child = tree.nodes[i].children[1]
            right_child = tree.nodes[i].children[2]
            left_prob = tree.nodes[i].rho
            right_prob = 1.0-tree.nodes[i].rho
    
            N_i = cur.num_leaves
            psi_term += (left_child.num_leaves-1) * log(left_prob)
            psi_term += (right_child.num_leaves-1) * log(right_prob)
            psi_term += -log(N_i) - log(N_i-1)

            N_i = cur.num_leaves
            p_s = 2/(N_i+1)

            rhot = tree.nodes[i].rhot
            if rhot == 1.0
                psi_term += log(p_s)
            else
                psi_term += log(1-p_s) + log(p_s) + (p_s-1)*log(rhot)
            end
            
        end
        if i != root.index
            parent_t = t[parent.index]
            inf_term += poisson_logpdf(tree.nodes[i].state, lambda * (parent_t - t[i]))
        end

    end

#    if inf_term < -3.0
#        println("INF_TERM: $inf_term")
#    end

    total_prob + psi_term + inf_term
end


function likelihood(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState,
                    indices::Array{Int64,1})
    YY = data.Ytrain
    (N,N) = size(YY[1])
    tree = model.tree
    Z = ConstructZ(tree)
    W = model.weights
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c
 
    total_prob = 0.0
    total_test_prob = 0.0

    for s = 1:length(data.Ytrain)
        Y = data.Ytrain[s]
        Ytest = data.Ytest[s]
        for i = indices
            for j = 1:N

    #            if i == 1 && mod(j, 20) == 0
    #                println("likelihood latent effect: ", squeeze(Z[i,:] * W * Z[j,:]')[1] )
    #            end
                oe = compute_observed_effects(model, model_spec, data, i, j)
                logit_arg = (Z[i,:] * W * Z[j,:]' + oe)[1]
                            
                
                total_prob += log_logit(logit_arg, Y[i,j])
                total_test_prob += log_logit(logit_arg, Ytest[i,j])
            end
        end
    end
    (total_prob, total_test_prob)
end

function test_likelihood_ij(model::ModelState,
                            model_spec::ModelSpecification,
                            data::DataState,
                            Z, #sparse or full binary matrix
                            i::Int,
                            j::Int)

    YYtest = data.Ytest
    (N,N) = size(YYtest[1])
    tree = model.tree
    W = model.weights
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c
 
    total_prob = 0.0
    total_test_prob = 0.0
  
    oe = compute_observed_effects(model, model_spec, data, i, j)
    logit_arg = (Z[i,:] * W * Z[j,:]' + oe)[1]
    LL = 0.0
    for s = 1:length(YYtest)
        Ytest = YYtest[s]
        LL += log_logit(logit_arg, Ytest[i,j])
    end
    (LL,logit_arg)
end

# Local pdfs for sampling

###################################
###### pdfs for psi updates #######
###################################

#@profile begin
#pdf for updating psi, the tree structure, assumes tree is already pruned at prune_index
function psi_infsites_logpdf(model::ModelState,
                             data::DataState,
                             pruned_index::Int64,
                             path::Array{Int64,1})
    (N,N) = size(data.Ytrain[1])
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
    pruned_indices = GetLeafToRootOrdering(tree, pruned_index)

    mu_1 = ones(2N - 1)
    mu_2 = ones(2N - 1)
    mu_3 = ones(2N - 1)

    pruned_node = tree.nodes[pruned_index]
    pruned_parent = pruned_node.parent

    self_direction = find(pruned_parent.children .== pruned_node)[1]

    pruned_mu_prop = self_direction == 1 ? (1-pruned_parent.rho) : pruned_parent.rho
    pruned_mu_prop *= pruned_parent.rhot

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    
        num_leaves = cur.num_leaves

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            cur_mu_prop *= parent.rhot
            mu_1[i] = mu_1[parent.index]*cur_mu_prop
        end

        mu_2[i] = mu_1[i] * pruned_mu_prop 

    end
    for j = reverse(pruned_indices)
        cur = tree.nodes[j]
        parent = cur.parent    
        if j == pruned_index
            mu_3[parent.index] = 1.0
        end

        self_direction = find(parent.children .== cur)[1]
        cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
        cur_mu_prop *= parent.rhot
        mu_3[j] = mu_3[parent.index]*cur_mu_prop
    end

    t_1 = mu_1.^gam
    t_2 = mu_2.^gam
    t_3 = mu_3.^gam

    t_1[1:N] = 0.0
    t_2[1:N] = 0.0
    t_3[1:N] = 0.0

    subtree_probs = -Inf * ones(2N - 1)
    for i = indices
        graft_node = tree.nodes[i]

        if i > N 
            t_graft = t_1[i]
        else
            parent = graft_node.parent
            self_direction = find(parent.children .== graft_node)[1]
            mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho 
            t_graft = t_1[parent.index] * (mu_prop * parent.rhot)^gam
        end

        if length(pruned_indices) > 0
            subtree_probs[i] = 0.0
        end
        for j = pruned_indices
            parent_t = t_3[tree.nodes[j].parent.index]
            cur = tree.nodes[j]
            poisson_mean = model.lambda * t_graft * (parent_t - t_3[j])
            subtree_probs[i] += poisson_logpdf(cur.state, poisson_mean)
        end
    end

    for i = indices
        cur = tree.nodes[i]

        if i == root.index
            parent_t1 = 1.0
        else
            parent_t1 = t_1[cur.parent.index]
        end

        parent_t2 = parent_t1 * pruned_mu_prop^gam
        poisson_mean_before = model.lambda * (parent_t1 - t_1[i])
        poisson_mean_after = model.lambda * (parent_t2 - t_2[i])
        descendant_mutation_probs[i] = poisson_logpdf(cur.state, poisson_mean_after) -
                                       poisson_logpdf(cur.state, poisson_mean_before)

        child_mutation_contributions = 0.0
        if i > N
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

        if i > N 
            pruned_parent_t = t_1[i]
        else
            parent = cur.parent
            self_direction = find(parent.children .== cur)[1]
            mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho 
            pruned_parent_t = parent_t1 * (mu_prop * parent.rhot)^gam
        end
        poisson_mean_after_above = model.lambda * (parent_t1 - pruned_parent_t)
        poisson_mean_after_below = model.lambda * (pruned_parent_t - t_2[i])

        if contains(path, i)
            (U, V) = all_splits([features_start:features_end])
            for j = 1:length(U)
                u = U[j]
                v = V[j]

                prob = poisson_logpdf(length(u), poisson_mean_after_below) +
                       poisson_logpdf(length(v), poisson_mean_after_above) -
                       poisson_logpdf(cur.state, poisson_mean_before)


                push!(total_probs, prob + child_mutation_contributions + psi_probs[i] + subtree_probs[i])
                if isnan(total_probs[end])
                    prob1 = poisson_logpdf(length(u), poisson_mean_after)
                    prob2 = poisson_logpdf(length(v), poisson_mean_before)
                    prob3 = -poisson_logpdf(cur.state, poisson_mean_before)
                    println("mean_after: $poisson_mean_after")
                    println("mean_before: $poisson_mean_before")
                    println("(u,v,i): $(length(u)), $(length(v)), $(cur.state)")
                    println("probs (u,v,i): $prob1, $prob2, $prob3")
                    println("prob ", prob)
                    println("child contributions ", child_mutation_contributions)
                    println("psi prob ", psi_probs[i])
                    println("subtree_probs ", subtree_probs[i])
                end
                # (node above which to attach, features moved below graft point, features moved above)
                push!(tree_states, (i, u, v) )                
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
                               data::DataState,
                               pruned_index::Int64,
                               path::Array{Int64,1})

    YY = data.Ytrain
    (N,N) = size(data.Ytrain[1])
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
            observed_parenthood_effects[i,l2] = compute_observed_effects(model, model_spec, data, l1, l2 )
            observed_childhood_effects[i,l2] = compute_observed_effects(model, model_spec, data, l2, l1)
        end
    end

    for i = subtree_indices
        cur = tree.nodes[i]
        features_start = weight_indices[i]
        features_end = features_start + cur.state - 1

        cur_features = [features_start:features_end] 

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

    for i = 1:length(nzI)
        push!(leaf_features[nzI[i]], nzJ[i])
    end

    latent_effects = zeros(Float64, (length(path), N, N))
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

            l1 = subtree_leaves[j1]

            if path_index == 1

                assert(cur.parent == Nil())
                append!(target_subtree_features[path_index,j1], subtree_features[j1])

                for j2 = 1:length(leaves)
                    l2 = leaves[j2]
                    assert(subtree_leaves[j1] != l2)
                    latent_effects[path_index,l1,l2] = sum(W[target_subtree_features[path_index,j1], leaf_features[l2]])
                    latent_effects[path_index,l2,l1] = sum(W[leaf_features[l2], target_subtree_features[path_index,j1]])
                end
            else
                features_start = weight_indices[parent.index]
                features_end = features_start + parent.state - 1

                new_features = [features_start:features_end]
                target_subtree_features[path_index,j1] = copy(target_subtree_features[parent_path_index,j1])
                append!(target_subtree_features[path_index,j1], new_features)

                for j2 = 1:length(leaves)
                    l2 = leaves[j2]
                    latent_effects[path_index,l1,l2] = sum(W[new_features, leaf_features[l2]]) +
                                                       latent_effects[parent_path_index,l1,l2]
                    latent_effects[path_index,l2,l1] = sum(W[leaf_features[l2], new_features]) +
                                                       latent_effects[parent_path_index,l2,l1]
                end
            end

        end

        for j1 = 1:length(subtree_leaves)
            l1 = subtree_leaves[j1]
            for j2 = 1:length(subtree_leaves)
                l2 = subtree_leaves[j2]
                latent_effects[path_index,l1,l2] = sum(W[target_subtree_features[path_index,j1], target_subtree_features[path_index,j2]])
            end
        end

        features_start = weight_indices[i]
        features_end = weight_indices[i] + cur.state - 1

        (V, T) = all_splits([features_start:features_end])

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
                    latent_effect_1 = sum(W[t, leaf_features[l2]]) + latent_effects[path_index,l1,l2]
                    latent_effect_2 = sum(W[leaf_features[l2], t]) + latent_effects[path_index,l2,l1]

#                    if l1 == 1 && mod(l2,20) == 0
#                        println("latent effect: ", latent_effect)
#                        println("leaf features: ", leaf_features[l2])
#                    end

                    for s = 1:length(YY)
                        Y = YY[s]
                        likelihood += log_logit(latent_effect_1 + observed_parenthood_effects[j1,l2], Y[l1,l2])
                        likelihood += log_logit(latent_effect_2 + observed_childhood_effects[j1,l2], Y[l2,l1])
                    end

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
                    latent_effect = sum(W[t, subtree_leaf_features2]) + sum(W[t, subtree_leaf_features1]) + sum(W[t,t]) + latent_effects[path_index,l1,l2]
                    for s = 1:length(YY)
                        Y = YY[s]
                        likelihood += log_logit(latent_effect + observed_parenthood_effects[j1,l2], Y[l1,l2])
                        #likelihood += log_logit(latent_effect + observed_childhood_effects[j1,l2], Y[l2,l1])
                    end

                    if likelihood == -Inf
                        println("subtree leaves")
                        println(latent_effect,",",observed_parenthood_effects[j1,l2],",",observed_childhood_effects[j1,l2])
                    end
                
                end
    
            end
            insert!(likelihoods, 1, likelihood)
            insert!(tree_states, 1, (i, rV[m], t))
        end       

    end

    (likelihoods, tree_states)
end


function prior_tree(tree::Tree,
                    pruned_index::Int)
    N::Int = (length(tree.nodes) + 1) / 2

    prob_diffs = -Inf*ones(2N - 1)
    prob_left = -Inf*ones(2N - 1)
    prob_right = -Inf*ones(2N - 1)
    rhot_prob = -Inf*ones(2N-1)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    i = 1
    while contains(subtree_indices, i)
        i += 1
    end
    root = FindRoot(tree, i)

    # indices will contain a leaves to root ordering of nodes in the pruned tree
    indices = GetLeafToRootOrdering(tree, root.index)

    pruned_indices = GetLeafToRootOrdering(tree, pruned_index)
    pruned_node = tree.nodes[pruned_index]
    pruned_parent = pruned_node.parent

    self_direction = find(pruned_parent.children .== pruned_node)[1]

    pruned_mu_prop = self_direction == 1 ? (1-pruned_parent.rho) : pruned_parent.rho

    pruned_leaves = pruned_node.num_leaves
    
    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    
        num_leaves = cur.num_leaves

        if i == root.index
            parent_rhot_prob = 0.0
            cur_prob = 0.0
        else
            parent_rhot_prob = rhot_prob[parent.index]
            self_direction = find(parent.children .== cur)[1]
            cur_prob = self_direction == 1 ? prob_left[parent.index] : prob_right[parent.index]
        end

        if i > N
            N_l = cur.children[1].num_leaves
            N_r = cur.children[2].num_leaves
            N_i = cur.num_leaves
            assert(N_l + N_r == N_i)

            ps_new = 2/(N_i+pruned_leaves+1)
            ps_old = 2/(N_i+1)

            rhot_prob[i] = cur.rhot == 1.0 ? log(ps_new)-log(ps_old) : 
                                (log(1-ps_new) + log(ps_new) + (ps_new-1)*log(cur.rhot) -
                                 log(1-ps_old) - log(ps_old) - (ps_old-1)*log(cur.rhot))
            rhot_prob[i] += parent_rhot_prob

            prob_right[i] = pruned_leaves * log(1-cur.rho) -
                            log(N_i+pruned_leaves) - log(N_i+pruned_leaves-1) +
                            log(N_i) + log(N_i-1) +
                            cur_prob
            prob_left[i] = pruned_leaves * log(cur.rho) - 
                           log(N_i+pruned_leaves) - log(N_i+pruned_leaves-1) +
                           log(N_i) + log(N_i-1) +
                           cur_prob
        end

        total_leaves = num_leaves + pruned_leaves
        p_s = 2/(total_leaves + 1)
        rhot = pruned_parent.rhot
        cur_rhot_prob = rhot == 1.0 ? log(p_s) : log(1-p_s) + log(p_s) + (p_s-1)*log(rhot) 

        prob_diffs[i] = cur_prob + parent_rhot_prob + cur_rhot_prob +
                        (num_leaves - 1) * log(pruned_mu_prop) -
                        log(total_leaves) - log(total_leaves-1) 
#        prob_above[i] = num_leaves * log(pruned_mu_prop)
#        prob_below[i] = pruned_leaves * log(mu_1[i])

#        if i != root.index
#            cum_prob_above[i] = cum_prob_above[parent.index] + prob_above[i]
#        end
    end


    prob_diffs
end


###################################
###### pdfs for rho, rhot updates #
###################################

function rhot_logf(rhot, l, r, N, gam, lambda, nu_l, nu_r,
                   k_l, k_r, K_l, K_r, T_l, T_r, S_l, S_r)
    Rl = l > N ? 1.0-(nu_l*rhot)^gam : 1.0
    Rr = r > N ? 1.0-(nu_r*rhot)^gam : 1.0
    k_l*log(Rl) + k_r*log(Rr) + gam*(K_l+K_r)*log(rhot) -
    lambda*( T_l*Rl + T_r*Rr + 
            (S_l*nu_l^gam + S_r*nu_r^gam)*rhot^gam )
end

function rhot_splits(rhot, p_s, l, r, N, gam, lambda, nu_l, nu_r, 
                     k_l, k_r, K_l, K_r, T_l, T_r, S_l, S_r)
    f1 = log(p_s) + 
         rhot_logf(1.0, l, r, N, gam, lambda, nu_l, nu_r, k_l, k_r, K_l, K_r, T_l, T_r, S_l, S_r)

    f2 = log(1-p_s)+log(p_s)+(p_s-1)*log(rhot) + 
         rhot_logf(rhot, l, r, N, gam, lambda, nu_l, nu_r, k_l, k_r, K_l, K_r, T_l, T_r, S_l, S_r)
   
    return [f1, f2] 
end


function rho_logpdf(rho, l, r, N, gam, lambda, nut_l, nut_r,
                    k_l, k_r, K_l, K_r, T_l, T_r, S_l, S_r, N_l, N_r)

    Rl = l > N ? 1.0-(nut_l*rho)^gam : 1.0
    Rr = r > N ? 1.0-(nut_r*(1-rho))^gam : 1.0

    k_l*log(Rl) + k_r*log(Rr) + 
    (gam*K_l+N_l)*log(rho) + (gam*K_r+N_r)*log(1-rho) -
    lambda*( T_l*Rl + T_r*Rr + 
             S_l*(nut_l*rho)^gam + S_r*(nut_r*(1-rho))^gam )
end


###################################
###### pdfs for W updates #########
###################################

# pdf for updating one element of W
function W_local_logpdf(model::ModelState,
                        data::DataState,
                        relevant_pairs::Array{Int64, 2},
                        latent_effects::Array{Float64, 2},
                        observed_effects::Array{Float64, 2},
                        w_old::Float64,
                        w_new::Float64)
    YY = data.Ytrain
    sigma = model.w_sigma 

    (_, npairs) = size(relevant_pairs)
  
    logprob = 0.0 
    for p = 1:npairs
        i = relevant_pairs[1,p]
        j = relevant_pairs[2,p]

        le = latent_effects[i,j]
        oe = observed_effects[i,j]

        le += w_new - w_old

        for s = 1:length(YY)
            Y = YY[s]   
            logprob += log_logit(le + oe, Y[i,j])
        end
    end 

    return logprob + normal_logpdf(w_new, sigma)
end

function log_logistic{T<:FloatingPoint}(effect::T, y::T)
    value = 0.0
    value = -(one(T)-y)*effect - log(one(T)+exp(-effect))
    if isinf(value)
        value = 0.0
    end
    return value
end

function log_logistic_dx(effect, y)
    -(1-y) + exp(-effect)/(1+exp(-effect))
end

#type LogLogistic <: BinaryFunctor end
#evaluate(::LogLogistic, x::Number, y::Number) = log_logistic(x,y)
#result_type(::LogLogistic, t1::Type, t2::Type) = promote_type(t1,t2)


# pdf for updating all elements of W
function W_full_logpdf(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       Z::Matrix{Int},
                       observed_effects::Array{Float64, 2})
    YY = data.Ytrain
    sigma = model.w_sigma

    if model_spec.positive_W
        W = exp(model.weights)
    else
        W = model.weights
    end

    latent_effects = Z*W*Z' 
    effects = latent_effects + observed_effects

    log_probs = 0.0 

    for p = 1:length(YY)
        Y = YY[p]

        log_probs_t = broadcast(log_logistic, effects, Y) 
        n_inds = find(Y .< 0)
        log_probs_t[n_inds] = 0.0
        log_probs += log_probs_t
    end

    if model_spec.positive_W # W is represented in log space 
        return sum(sum(log_probs)) + sum(exp_logpdf(W, sigma)) + sum(model.weights)
    else
        return sum(sum(log_probs)) + sum(normal_logpdf(W, sigma))
    end
end

function W_full_logpdf_gradient(model::ModelState,
                                model_spec::ModelSpecification,
                                data::DataState,
                                Z::Matrix{Int},
                                observed_effects::Array{Float64, 2})

    YY = data.Ytrain
    sigma = model.w_sigma

    if model_spec.positive_W
        W = exp(model.weights)
    else
        W = model.weights
    end

     
    latent_effects = Z*W*Z' 

    effects = latent_effects + observed_effects

    log_sigmoid_gradient = 0.0
    for p = 1:length(YY)
        Y = YY[p]

        gradient = broadcast(log_logistic_dx, effects, Y)
        n_inds = find(Y .< 0)
        gradient[n_inds] = 0.0
        log_sigmoid_gradient += gradient
    end

    ZLZ = Z'*log_sigmoid_gradient*Z

    if model_spec.symmetric_W
        ZLZ_L = tril(ZLZ)
        ZLZ_L[diagind(ZLZ)] = 0.0
        ZLZ += ZLZ_L'
        symmetrize!(ZLZ)
    end

    if model_spec.positive_W # W is represented in log space 
        LL_grad = ZLZ.*W
        return LL_grad + exp_logpdf_dx(W, sigma) + ones(size(W))
    else
        LL_grad = ZLZ.*W
        return LL_grad + normal_logpdf_dx(W, sigma)
    end

end
###################################
###### pdfs for intercept updates #
###################################


function intercept_logpdf(model::ModelState,
                     model_spec::ModelSpecification,
                     data::DataState,
                     latent_effects::Array{Float64, 2},
                     observed_effects::Array{Float64, 2},
                     old_intercept::Float64,
                     new_intercept::Float64)

    YY = data.Ytrain
    (N,N) = size(YY[1])
     
    b_sigma = model.b_sigma
    total_prob = 0.0
    total_prob += normal_logpdf(new_intercept, b_sigma)

    for i = 1:N
        for j = 1:N
            oe = observed_effects[i,j] + new_intercept - old_intercept
            logit_arg = latent_effects[i,j] + oe
            for s = 1:length(YY)
                Y = YY[s]
                total_prob += log_logit(logit_arg, Y[i,j])
            end
        end
    end
    total_prob
end

function ab_logpdf(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState,
                  latent_effects::Array{Float64, 2})

    YY = data.Ytrain
    (N,N) = size(YY[1])
     
    b_sigma = model.b_sigma
    a = model.a
    b = model.b

    observed_effects = construct_observed_effects(model, model_spec, data)
    effects = latent_effects + observed_effects


    total_prob = 0.0
    a_gradient = 0.0
    b_gradient = 0.0

    log_sigmoid_gradient = 0.0
    for p = 1:length(YY)
        Y = YY[p]
        gradient = broadcast(log_logistic_dx, effects, Y)
        n_inds = find(Y .< 0)
        gradient[n_inds] = 0.0
        log_sigmoid_gradient += gradient

        total_prob += sum(broadcast(log_logit, effects, Y))
    end

    if model_spec.use_parenthood
        total_prob += sum(normal_logpdf(a, b_sigma))
        a_gradient += normal_logpdf_dx(a, b_sigma)
        a_gradient += squeeze(sum(log_sigmoid_gradient,1),1)
    end

    if model_spec.use_childhood
        total_prob += sum(normal_logpdf(b, b_sigma))
        b_gradient += normal_logpdf_dx(b, b_sigma)
        b_gradient += squeeze(sum(log_sigmoid_gradient,2),2)
    end


    total_prob, a_gradient, b_gradient
end
###################################
###### pdfs for vardim updates ####
###################################
#@profile begin


function vardim_local_logpdf(model::ModelState,
                             data::DataState,
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
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            # W_local_logpdf includes the prior term for w_old as if it were normal,
            # need to subtract it out here
            logprob += aug_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    elseif u == model.tree.nodes[node_index].state
        if w_is_auxiliary
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            logprob += aug_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    else
        logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                  observed_effects, w_old, w_new)
    end

    logprob + poisson_logpdf(u, effective_lambda)
end

function vardim_local_splits(model::ModelState,
                             model_spec::ModelSpecification,
                             data::DataState,
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
        logprob = vardim_local_logpdf(model, data, relevant_pairs, latent_effects,
                                      observed_effects, mixture_component, 
                                      u[mixture_component], node_index,
                                      effective_lambda, w_old, w_new,
                                      w_is_auxiliary[mixture_component]) 
        logprobs[mixture_component] = logprob

    end
    logprobs += vardim_multiplier_terms(model_spec, u)
    logprobs + constant_terms
end

function vardim_local_sum(model::ModelState,
                          model_spec::ModelSpecification,
                          data::DataState,
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
    logprobs = vardim_local_splits(model, model_spec, data, relevant_pairs, latent_effects,
                                   observed_effects, u, weight_index_pointers,
                                   node_index, effective_lambda, constant_terms, 
                                   w_old, w_new, w_is_auxiliary)
    logsumexp(logprobs)
end



function vardim_logpdf(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       observed_effects::Array{Float64, 2},
                       u::Int64,
                       L::Int64,
                       aug_k::Int64,
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64)


    tmodel = copy(model)
    tmodel.tree.nodes[node_index].state = L+1
    Z = ConstructZ(tmodel.tree)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + L - 1
    new_index = end_index + 1 

    if model_spec.positive_W
        W = exp(model.weights)
    else
        W = copy(model.weights)
    end

    if u < L
        assert( aug_k < new_index)
        W[:,new_index] = 0.0
        W[new_index,:] = 0.0
        W[aug_k,:] = 0.0
        W[:,aug_k] = 0.0
    elseif L == u
        W[:,new_index] = 0.0
        W[new_index,:] = 0.0
    elseif u > L
        assert(aug_k == new_index)
    end

    if model_spec.symmetric_W
        W = copy(W)
        symmetrize!(W)
    end

    YY = data.Ytrain
    (N,N) = size(YY[1])

    logprob = 0.0

    gam = model.gamma
    lambda = model.lambda
    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    ZWZ = Z*W*Z' 
    effects = ZWZ + observed_effects

    logprob = 0.0
    gradient = 0.0
    for s = 1:length(YY)
        Y = YY[s]
        LL = broadcast(log_logit, effects, Y)
        LL_grad = broadcast(log_logistic_dx, effects, Y)
        missing_data = find(Y .< 0)
        LL[missing_data] = 0.0
        LL_grad[missing_data] = 0.0
      
        logprob += sum(LL) 
        gradient += LL_grad 
    end

    ZLZ = Z'*gradient*Z

    if model_spec.symmetric_W
        ZLZ_L = tril(ZLZ)
        ZLZ_L[diagind(ZLZ)] = 0.0
        ZLZ += ZLZ_L'
        symmetrize!(ZLZ)
    end


    gradient = ZLZ.*W

    W_logpdf = model_spec.W_logpdf
    W_logpdf_gradient = model_spec.W_logpdf_gradient
    W_full = model.weights
    (K,K) = size(W)
    sigma = model.w_sigma

    W_logprob = W_logpdf(W, sigma)    

    if model_spec.diagonal_W
        logprob += sum(diag(W_logprob)) 
    elseif model_spec.symmetric_W
        logprob += sum(triu(W_logprob))
    else
        logprob += sum(W_logprob)
    end

    W_grad = W_logpdf_gradient(W, sigma)
    gradient += W_grad

    logprob += poisson_logpdf(u, effective_lambda)

    return logprob, gradient 
end

function vardim_splits(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       observed_effects::Array{Float64, 2},
                       U::Array{Int64, 1},
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64)
    logprobs = zeros(length(U))
    gradients = cell(length(U))
    assert(length(U) > 1)
    L = U[end-1] 

    start_index = weight_index_pointers[node_index]
    end_index = start_index + L - 1
    new_index = end_index + 1 

    for u_ind = 1:length(U)
        u = U[u_ind]
          
        if u < L 
            aug_k = start_index + u_ind - 1
        elseif u == L
            aug_k = 0
        elseif u > L
            aug_k = new_index
        end

        logprobs[u_ind], gradients[u_ind] = 
                              vardim_logpdf(model, model_spec, data,
                              observed_effects, u, L, aug_k, 
                              weight_index_pointers, node_index,
                              effective_lambda)

    end
    logprobs += vardim_multiplier_terms(model_spec, U)
    logprobs, gradients
end

function vardim_sum(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState,
                    observed_effects::Array{Float64, 2},
                    U::Array{Int64, 1},
                    weight_index_pointers::Array{Int64, 1},
                    node_index::Int64,
                    effective_lambda::Float64)

    logprobs, gradients = vardim_splits(model, model_spec, data, observed_effects,
                                        U, weight_index_pointers, node_index,
                                        effective_lambda)

    logprob = logsumexp(logprobs)
    gradient = logsumexp_d_dx(logprobs, gradients)

    logprob, gradient
end

function vardim_multiplier_terms(model_spec::ModelSpecification,
                                 u::Array{Int64,1})
    
    if length(u) == 1 # the case L < 0
        return [log(model_spec.rrj_jump_probabilities[1])]
    end

    logprobs = zeros(length(u))
    logprobs += log(model_spec.rrj_jump_probabilities[u - u[1] + 1])
    # for u with K < L == u[end-1], we need the 1/L term
    logprobs -= [u[x] < u[end-1] ? log(u[end-1]) : 0.0 for x = 1:length(u)]

    logprobs
end

#end #profile
# Utility Functions

# Assumes model.augmented_weights is up to date!
function compute_constant_terms(model::ModelState,
                                model_spec::ModelSpecification,
                                data::DataState,
                                relevant_pairs::Array{Array{Int64, 2}, 2},
                                latent_effects::Array{Array{Float64, 2}, 1},
                                observed_effects::Array{Float64, 2},
                                U::Array{Int64, 1},
                                weight_index_pointers::Array{Int64, 1},
                                node_index::Int64,
                                effective_lambda::Float64,
                                total_logprob::Float64)

    assert(length(U) > 0)
    logprobs = total_logprob * ones(length(U))

    current_ind = length(U) - 1

    current_u = U[current_ind]
    W = model.weights

    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1


    for u_ind = 1:length(U)
        u = U[u_ind]

        aug_k = u < current_u ? start_index + u_ind - 1 : u > current_u ? new_index : 0

        prob_tuple = adjust_model_logprob(model, model_spec, data, relevant_pairs,
                         latent_effects[current_ind], latent_effects[u_ind],
                         observed_effects, aug_k, current_u, u, weight_index_pointers,
                         node_index, effective_lambda)


        logprobs[u_ind] += sum(prob_tuple)
#        if model_spec.debug
#            println("u_ind, prob_tuple: ", (u_ind, prob_tuple))
#            println("u_ind, u, current_u, aug_k: ", (u_ind, u, current_u, aug_k))
#        end
    end 
    logprobs
end

# Adjust model_logprob when L != u
# Assumes model.augmented_weights is up to date!
function adjust_model_logprob(model::ModelState,
                              model_spec::ModelSpecification,
                              data::DataState,
                              relevant_pairs::Array{Array{Int64, 2}, 2},
                              current_latent_effects::Array{Float64, 2},
                              latent_effects::Array{Float64, 2},
                              observed_effects::Array{Float64, 2},
                              aug_k::Int,
                              u::Int,
                              L::Int,
                              weight_index_pointers::Array{Int64, 1},
                              node_index::Int64,
                              effective_lambda::Float64)

    if L < 0
        return (0.0, 0.0, 0.0)
    end

    YY = data.Ytrain

    augW = model.augmented_weights
    aug_inds, _ = get_augmented_submatrix_indices(augW, node_index, 1)

    W = augW[aug_inds, aug_inds]
    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1

    aug_terms = 0.0
    prior_terms = 0.0
    likelihood_terms = 0.0
    if model_spec.diagonal_W
        if L < u
            removed_ind = aug_k
            aug_terms += aug_logpdf(W[new_index,new_index], model.nu)
            aug_terms += aug_logpdf(W[removed_ind, removed_ind], model.nu)
            prior_terms -= normal_logpdf(W[removed_ind, removed_ind], model.w_sigma)

            rpairs = relevant_pairs[removed_ind, removed_ind]
            for p = 1:size(rpairs)[2]
                i = rpairs[1,p]
                j = rpairs[2,p]

                le = latent_effects[i,j]
                oe = observed_effects[i,j]
                old_le = current_latent_effects[i,j]
                for s = 1:length(YY)
                    Y = YY[s]
                    likelihood_terms += log_logit(le + oe, Y[i,j])
                    likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                end
            end 
        elseif L == u
            aug_terms += aug_logpdf(W[new_index,new_index], model.nu)
        else 
            prior_terms += normal_logpdf(W[new_index,new_index], model.w_sigma)
            rpairs = relevant_pairs[new_index, new_index]
            for p = 1:size(rpairs)[2]
                i = rpairs[1,p]
                j = rpairs[2,p]

                le = latent_effects[i,j]
                oe = observed_effects[i,j]
                old_le = current_latent_effects[i,j]
                for s = 1:length(YY)
                    Y = YY[s]
                    likelihood_terms += log_logit(le + oe, Y[i,j])
                    likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                end 
            end 
        end
    else
        for k1 = 1:K # start_index:new_index
            k2_range = k1 != new_index ? new_index :
                       model_spec.symmetric_W ? (1:k1) : 
                       (1:K)
            for k2 = k2_range
                if L <= u
                    aug_terms += aug_logpdf(W[k1,k2], model.nu)
                else
                    prior_terms += normal_logpdf(W[k1,k2], model.w_sigma)
                end
            end
        end

        c_latent_effects = copy(current_latent_effects)

        if L < u || L > u
            assert( L < u || aug_k == new_index)
            for k1 = 1:K
                if model_spec.symmetric_W
                    k2_range = k1 == aug_k ? (1:k1) :
                               k1 < aug_k ? aug_k : []
                else
                    k2_range = k1 == aug_k ? (1:K) : aug_k
                end
                for k2 = k2_range
                    if k1 != new_index && k2 != new_index
                        aug_terms += aug_logpdf(W[k1,k2], model.nu)
                        prior_terms -= normal_logpdf(W[k1,k2], model.w_sigma)
                    end

                    rpairs = relevant_pairs[k1,k2]
                    for p = 1:size(rpairs)[2]
                        i = rpairs[1,p]
                        j = rpairs[2,p]

                        le = latent_effects[i,j]
                        oe = observed_effects[i,j]
                        old_le = c_latent_effects[i,j]
                        # dirty hack to prevent double counting (ie from visiting the same
                        # pair more than once)
                        c_latent_effects[i,j] = le
                        for s = 1:length(YY)
                            Y = YY[s]
                            likelihood_terms += log_logit(le + oe, Y[i,j])
                            likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                        end
                    end
                end
            end
        end
    end

    prior_terms += poisson_logpdf(L,effective_lambda) - 
                   poisson_logpdf(u,effective_lambda)

    (prior_terms, likelihood_terms, aug_terms)
end



function adjust_component_latent_effects(model::ModelState,
                                         model_spec::ModelSpecification,
                                         latent_effects::Array{Float64,2},
                                         relevant_pairs::Array{Array{Int64,2},2},
                                         u::Int,
                                         start_index::Int,
                                         end_index::Int)

    if u < 0
        component_latent_effects = Array(Array{Float64,2}, 1)
        component_latent_effects[1] = latent_effects
        return component_latent_effects
    end

    W = model.weights
    (K,K) = size(W)
    (rpK,rpK) = size(relevant_pairs)
#    println("rpK: ", rpK)
#    println("K: ", K)
    assert( rpK == K )

    new_k = end_index + 1
    component_latent_effects = Array(Array{Float64,2}, u+2)


    for u_ind = 1:u+2
        if u_ind <= u
            removed_k = start_index + u_ind - 1
            assert(removed_k <= K)

            component_latent_effects[u_ind] =
                adjust_latent_effects(model, model_spec, latent_effects, relevant_pairs,
                                      u, u-1, removed_k,new_k)
        elseif u_ind == u+1
            component_latent_effects[u_ind] = copy(latent_effects)
        elseif u_ind == u+2
            assert(new_k <= K)
            component_latent_effects[u_ind] =
                adjust_latent_effects(model, model_spec, latent_effects, relevant_pairs,
                                      u, u+1, new_k, new_k)
        end
    end

    component_latent_effects
end

function compute_component_latent_effects(model::ModelState,
                                          model_spec::ModelSpecification,
                                          u::Int,
                                          start_index::Int,
                                          end_index::Int,
                                          node_index::Int)

    if u < 0
        component_latent_effects = Array(Array{Float64,2}, 1)
        component_latent_effects[1] = latent_effects
        return component_latent_effects
    end

    W = model.weights
    (K,K) = size(W)

    new_k = end_index + 1
    component_latent_effects = Array(Array{Float64,2}, u+2)


    for u_ind = 1:u+2
        if u_ind <= u
            removed_k = start_index + u_ind - 1
            assert(removed_k <= K)

            component_latent_effects[u_ind] =
                compute_latent_effects(model, model_spec, 
                                      u, u-1, removed_k, node_index)
        elseif u_ind == u+1
            component_latent_effects[u_ind] = 
                compute_latent_effects(model, model_spec, 
                                       u, u, 0, node_index)
        elseif u_ind == u+2
            assert(new_k <= K)
            component_latent_effects[u_ind] =
                compute_latent_effects(model, model_spec, 
                                      u, u+1, new_k, node_index)
        end
    end

    component_latent_effects
end


# adjusts latent effects when L != u 
function adjust_latent_effects(model::ModelState,
                               model_spec::ModelSpecification,
                               latent_effects::Array{Float64,2},
                               relevant_pairs::Array{Array{Int64,2},2},
                               u::Int,
                               L::Int,
                               aug_k::Int,
                               new_k::Int)

    W = model.weights
    (K,K) = size(relevant_pairs)

    assert(aug_k <= K)
    assert(new_k <= K)
    new_latent_effects = copy(latent_effects)
    if L == u - 1
        assert(aug_k < new_k || new_k == 0)

        removed_k = aug_k
        k_range0 = [1:K]
        # latent effects assumes new dimensions haven't been added, so don't remove them!
        if new_k != 0 # only need to worry about new dimension if we have one
            splice!(k_range0, new_k)
        end

        k_range = model_spec.diagonal_W ? removed_k : k_range0

        for k = k_range
            for p = 1:size(relevant_pairs[k, removed_k])[2]
                i = relevant_pairs[k, removed_k][1,p]
                j = relevant_pairs[k, removed_k][2,p]
                new_latent_effects[i,j] -= W[k, removed_k]

            end
            if k != removed_k
                for p = 1:size(relevant_pairs[removed_k, k])[2]
                    i = relevant_pairs[removed_k, k][1,p]
                    j = relevant_pairs[removed_k, k][2,p]
                    new_latent_effects[i,j] -= W[removed_k, k]
                end
            end
        end

    elseif L == u + 1
        new_k = aug_k
        k_range = model_spec.diagonal_W ? new_k : (1:K)
        
        for k = k_range

            if model_spec.symmetric_W
                kmin = min(k,new_k)
                kmax = max(k,new_k)

                if kmax > kmin
                    assert(size(relevant_pairs[kmin,kmax])[2] == 0)
                end
            end

            for p = 1:size(relevant_pairs[k, new_k])[2]
                i = relevant_pairs[k, new_k][1,p]
                j = relevant_pairs[k, new_k][2,p]
                new_latent_effects[i,j] += W[k, new_k]
            end
            if k != new_k
                for p = 1:size(relevant_pairs[new_k, k])[2]
                    i = relevant_pairs[new_k, k][1,p]
                    j = relevant_pairs[new_k, k][2,p]
                    new_latent_effects[i,j] += W[new_k, k]
                end
            end
        end

    end

    new_latent_effects
end

# compute latent effects from scratch
function compute_latent_effects(model::ModelState,
                                model_spec::ModelSpecification,
                                u::Int,
                                L::Int,
                                aug_k::Int,
                                node_index::Int;
                                include_gradient::Bool=false)

    tmodel = copy(model)
    tmodel.tree.nodes[node_index].state = L
    Z = ConstructZ(tmodel.tree)
    w_index_pointers = weight_index_pointers(model.tree)
    start_index = w_index_pointers[node_index]
    new_index = start_index + u 
    if L < u
        assert( aug_k < new_index)
        W = delete_row_and_col(model.weights, new_index)
        W = delete_row_and_col(W, aug_k)
    elseif L == u
        W = delete_row_and_col(model.weights, new_index)
    elseif L > u
        assert(aug_k == new_index)
        W = model.weights
    end

    if model_spec.symmetric_W
        W = copy(W)
        symmetrize!(W)
    end

    Z*W*Z'
end

function compute_all_relevant_pairs(model_spec::ModelSpecification,
                                    data::DataState,
                                    K::Int,
                                    Z)

    YY = data.Ytrain
    Ynn = zeros(size(YY[1]))
    Ynn[find(YY[1] .>= 0)] = 1
    Ynn = sparse(Ynn)
    relevant_pairs = [zeros(Int64,(0,0)) for x = 1:K, y = 1:K]

    for k1 = 1:K
        k2_range = model_spec.diagonal_W ? k1 :
                   model_spec.symmetric_W ? (1:k1) :
                   (1:K)
        for k2 = k2_range
            relevant_pairs[k1,k2] = compute_relevant_pairs(model_spec,Ynn,Z,k1,k2)
        end
    end 

    return relevant_pairs
end

# Compute relevant pairs for newly introduced weight parameters
# should be called without changing model.tree from the original
function compute_new_relevant_pairs(model::ModelState,
                                    model_spec::ModelSpecification,
                                    data::DataState,
                                    relevant_pairs::Array{Array{Int64,2},2},
                                    node_index::Int,
                                    start_index::Int,
                                    end_index::Int)
    (K,K) = size(relevant_pairs)
    K += 1
    N::Int = (length(model.tree.nodes)+1)/2
    new_k = end_index + 1

#    println("new_k: ", new_k)
#    println("K: ", K)
    assert(new_k <= K)

    tree = model.tree
    Z = ConstructZ(tree)
    # construct new Z col
    #println("construct Znew")
    leaves = GetLeaves(tree, node_index)
    Zk = zeros(Int64, N)
    Zk[leaves] = 1

    Znew = zeros(Int64, (N, K))
    Znew[:,1:end_index] = Z[:,1:end_index]
    Znew[:,end_index+1] = Zk
    Znew[:,end_index+2:end] = Z[:,end_index+1:end]        
    Znew = sparse(Znew)

    new_relevant_pairs = [zeros(Int64,(0,0)) for x = 1:K, y = 1:K]
    nonzero_element_indices = [x <= end_index ? x : x + 1 for x in 1:K-1]

    YY = data.Ytrain
    Ynn = zeros(size(YY[1]))
    Ynn[find(YY[1] .>= 0)] = 1
    Ynn = sparse(Ynn)

    new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = copy(relevant_pairs)
    k1_range = model_spec.diagonal_W ? (start_index:end_index+1) : (1:K)
    for k1 = k1_range
        k2 = new_k
        if model_spec.symmetric_W
            k11 = max(k1,k2)
            k22 = min(k1,k2)
            new_relevant_pairs[k11,k22] = compute_relevant_pairs(model_spec,Ynn,Znew,k11,k22)
        else
            new_relevant_pairs[k1,k2] = compute_relevant_pairs(model_spec,Ynn,Znew,k1,k2)
            new_relevant_pairs[k2,k1] = compute_relevant_pairs(model_spec,Ynn,Znew,k2,k1)
        end
    end
    new_relevant_pairs
end

function compute_unaugmented_prob(model::ModelState,
                                  model_spec::ModelSpecification,
                                  relevant_pairs::Array{Array{Int64, 2}, 2},
                                  latent_effects::Array{Array{Float64, 2}, 1},
                                  observed_effects::Array{Float64, 2},
                                  U::Array{Int64, 1},
                                  weight_index_pointers::Array{Int64, 1},
                                  node_index::Int64,
                                  effective_lambda::Float64,
                                  u_index::Int64,
                                  augmented_logprob::Float64)
    W = model.weights
    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1

    u = U[u_index]

    current_u = U[end-1]

#    println(U)
#    println(u_index)
    assert(current_u == model.tree.nodes[node_index].state)

    unaugmented_logprob = augmented_logprob
    if model_spec.diagonal_W
        if u < current_u
            removed_ind = start_index + u_index - 1
            unaugmented_logprob -= aug_logpdf(W[removed_ind,removed_ind], model.nu)
            unaugmented_logprob -= aug_logpdf(W[new_index,new_index], model.nu)
        elseif u == current_u 
            unaugmented_logprob -= aug_logpdf(W[new_index,new_index], model.nu)
        end #if u => current_u, no changes are necessary
    else
        if u < current_u
            removed_ind = start_index + u_index - 1
            for k1 = 1:K
                if model_spec.symmetric_W
                    @assert removed_ind < new_index
                    k2_range = k1 == removed_ind || k1 == new_index ? (1:k1) : 
                               k1 < removed_ind ? [] :
                               k1 < new_index ? removed_ind :
                               [removed_ind, new_index]
                    
                else
                    k2_range = k1 == removed_ind || k1 == new_index ? 
                               (1:K) : [removed_ind, new_index]
                end

                for k2 = k2_range
                    unaugmented_logprob -= aug_logpdf(W[k1,k2],model.nu)
                end
            end
        elseif u == current_u
            for k1 = 1:K
                if model_spec.symmetric_W
                    k2_range = k1 == new_index ? (1:k1) : 
                               k1 < new_index ? [] :
                               new_index
                else
                    k2_range = k1 == new_index ? (1:K) : new_index
                end
                for k2 = k2_range
                    unaugmented_logprob -= aug_logpdf(W[k1,k2],model.nu)
                end
            end
        end
    end 

    unaugmented_logprob
end

# Find all pairs i,j such that Z[i,k1]*Z[j,k2] = 1
function compute_relevant_pairs(model_spec::ModelSpecification,
                                mask, #sparse or full binary array, so we only inlude training data
                                Z, #sparse or full binary array
                                k1::Int64,
                                k2::Int64)
    ZZ = Z[:,k1]*Z[:,k2]'
    if model_spec.symmetric_W
        ZZ = ZZ + ZZ' - ZZ .* ZZ'
    end
    #ZZ = ZZ.*mask #only need pairs from Ytrain
    (I, J) = findn(ZZ)
    [I', J']
end

function symmetrize!(A::Array)
    triu_inds = find(triu(ones(size(A)),1))
    A[triu_inds] = tril(A,1)'[triu_inds];
end

# Doesn't quite handle the symmetric case yet (eg the a[i] and b[j])
function compute_observed_effects(model::ModelState,
                                  model_spec::ModelSpecification,
                                  data::DataState,
                                  i::Int64,
                                  j::Int64)

    X_r = data.X_r
    X_p = data.X_p
    X_c = data.X_c
    observed_effect = model.c

    if model_spec.use_pairwise
        observed_effect += model.beta' * X_r[i,j,:]
    end

    if model_spec.use_parenthood
        if length( X_p[j,:]) > 0
            observed_effect += model.beta_p' * X_p[j,:]
        end
        observed_effect += model.a[j]
    end

    if model_spec.use_childhood
        if length( X_c[i,:]) > 0
            observed_effect += model.beta_c' * X_c[i,:]
        end
        if model_spec.symmetric_W
            observed_effect += model.a[i]
        else
            observed_effect += model.b[i]
        end
    end                                  

    observed_effect
end

function construct_observed_effects(model::ModelState,
                                    model_spec::ModelSpecification,
                                    data::DataState)
    N = size(data.Ytrain[1],1) 
    observed_effects = zeros(Float64, (N,N))
    for i = 1:N
        for j = 1:N
            observed_effects[i,j] = compute_observed_effects(model, model_spec, data, i, j)
        end
    end
    observed_effects
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
                push!(u, v[j])
            end
            x = x >> 1
        end
        push!(U, u)
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

