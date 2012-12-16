load("model.jl")
load("tree.jl")
load("probability_util.jl")
load("pdf.jl")
load("slicesampler.jl")

function mcmc(Y::Array{Int64,2},
              X_r::Array{Float64,3},
              X_p::Array{Float64,2},
              X_c::Array{Float64,2},
              lambda::Float64,
              gamma::Float64,
              model_spec::ModelSpecification,
              iterations::Int64)

    (N,N) = size(Y)

    if model_spec.use_pairwise
        (N1, N2, pairwise_p) = size(X_r)
        assert(N1 == N2)
        assert(N == N2)
    end

    if model_spec.use_parenthood
        (N1, parent_p) = size(X_p)
        assert(N1 == N)
    end
    if model_spec.use_childhood
        (N2, children_p) = size(X_c)
        assert(N2 == N)
    end

    U = zeros(Int64,2N-1)
    for i = 1:2N-1
        u = rand()
        if u > .995
            U[i] = 2
        elseif u > .95
            U[i] = 1
        end
    end

    W = randn(sum(U),sum(U))
    tree = Tree(U)
    model = ModelState(lambda,gamma,1.0,1.0,tree,W,[0.0],[0.0],[0.0],[0.0],[0.0],0.0)


    #perform some unit tests
  
#    W_copy = copy(W) 
#    tree_copy = copy(tree)
#    assert(tree_copy == tree)
#
#    weight_indices = weight_index_pointers(tree)
#    prune_index = 1
#    for i = 1:length(tree.nodes)
#        if tree.nodes[i].parent.state > 0
#            prune_index = i
#            parent = tree.nodes[prune_index].parent
#            if parent.children[1].index == prune_index
#                sibling = parent.children[2]
#            else
#                sibling = parent.children[1]
#            end
#            if sibling != Nil() && sibling.state > 0
#                break
#            end
#        end
#    end
#    println("Testing prune: ", prune_index)
#    parent = tree.nodes[prune_index].parent
#    if parent.children[1].index == prune_index
#        sibling = parent.children[2]
#    else
#        sibling = parent.children[1]
#    end
#
#
#    st1 = sibling.state 
#
#    prune_tree!(model,prune_index)
#
#    st2 = sibling.state
#    println("st2 - st1: ", st2 - st1)
#    println("st1, st2: ", st1, ", ", st2)
#    sib_ind = sibling.index
#    weight_indices2 = weight_index_pointers(tree)
#
#    start_index = weight_indices2[sib_ind]
#    end_index = start_index + st2 - 1
#
#    parent_features = linspace(start_index + st1,end_index,st2 - st1)
#    sibling_features = linspace(start_index,start_index + st1 - 1,st1)
#    graft_tree!(model,prune_index, sib_ind, parent_features, sibling_features)
#
#    assert(tree_copy == tree)
#    assert(norm(W) == norm(model.weights))
#    assert(W == model.weights)

    ###

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prior = prior(model) 
        tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
        println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
        if iter == 10
            model_spec.debug = true
        end
        mcmc_sweep(model, model_spec, Y, X_r, X_p, X_c)
    end
end


function mcmc_sweep(model::ModelState,
                    model_spec::ModelSpecification,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2})

    (N,N) = size(Y)
    tree = model.tree
    Z = ConstructZ(tree)
    W = model.weights

    latent_effects = zeros(Float64, (N,N))
    observed_effects = zeros(Float64, (N,N))
    for i = 1:N
        for j = 1:N
            latent_effects[i,j] = (Z[i,:] * W * Z[j,:]')[1]
            observed_effects[i,j] = compute_observed_effects(model, model_spec, i, j, X_r, X_p, X_c)
        end
    end


    W_time = time()
    sample_W(model, model_spec, Y, X_r, X_p, X_c, Z, latent_effects, observed_effects)
    W_time = time() - W_time

    tree_prior = prior(model) 
    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)


    psi_time = time()
    sample_psi(model, model_spec, Y, X_r, X_p, X_c, latent_effects, observed_effects)
    psi_time = time() - psi_time

    tree_prior = prior(model) 
    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    vardim_time = time()
    sample_Z(model, model_spec, Y, X_r, X_p, X_c, latent_effects, observed_effects)
    vardim_time = time() - vardim_time



    (K,K) = size(model.weights)
    println("MCMC Timings (psi, W, vardim) = ", (psi_time, W_time, vardim_time))
    println("K: ", K)
end

function sample_W(model::ModelState,
                  model_spec::ModelSpecification,
                  Y::Array{Int64,2},
                  X_r::Array{Float64,3},
                  X_p::Array{Float64,2},
                  X_c::Array{Float64,2},
                  Z,#sparse or full 2d array
                  latent_effects::Array{Float64,2},
                  observed_effects::Array{Float64,2})

    (K,K) = size(model.weights)
    (N,N) = size(Y)
    W = model.weights

    num_W_sweeps = 5
    num_W_slice_steps = 1

   
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c

#    tree_prior = prior(model) 
#    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    for iter = 1:num_W_sweeps
   
        for k1 = 1:K
            println("W sampling k1: ", k1)
            for k2 = 1:K
                relevant_pairs = compute_relevant_pairs(Z,k1,k2)
                w_cur = W[k1,k2]
                w_old = w_cur
                g = x -> W_local_logpdf(model, Y, relevant_pairs, latent_effects, observed_effects, w_old, x)

                gx0 = g(w_cur)
                for slice_iter = 1:num_W_slice_steps
                    (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                end

                 
                W[k1,k2] = w_cur

                
                for p = 1:size(relevant_pairs)[2]
                    i = relevant_pairs[1,p]
                    j = relevant_pairs[2,p]

                    latent_effects[i,j] += w_cur - w_old
                end
#                if k1 == 45 && k2 == 40
#                    println(relevant_pairs)
#                    assert(false)
#                end

#                if model_spec.debug
#                    tree_prior = prior(model) 
#                    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#                    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
#
#                    for p = 1:size(relevant_pairs)[2]
#                        i = relevant_pairs[1,p]
#                        j = relevant_pairs[2,p]
#                        
#                        le = (Z[i,:] * W * Z[j,:]')[1]
#
#                        if abs(le - latent_effects[i,j]) > 10.0^-6
#                            println("p,i,j: ", (p,i,j))
#                            println("le, latent_effects: ", (le, latent_effects[i,j]))
#                            println("pairs: ", relevant_pairs)
#                            (I,J) = findn(Z[:,k1]*Z[:,k2]')
#                            println("I,J: ", (I,J))
#                            
#                            assert(le == latent_effects[i,j])
#                        end
#                    end
#
#                end


            end

        end

 
    end

end

function sample_Z(model::ModelState,
                  model_spec::ModelSpecification,
                  Y::Array{Int64,2},
                  X_r::Array{Float64,3},
                  X_p::Array{Float64,2},
                  X_c::Array{Float64,2},
                  latent_effects::Array{Float64,2},
                  observed_effects::Array{Float64,2})

    (N,N) = size(Y)
    tree = model.tree
    num_W_sweeps = 3
    num_W_slice_steps = 1
    # Sample new features from root to leaf
    root = FindRoot(tree, 1) 
    leaf_to_root = GetLeafToRootOrdering(tree, root.index)
    root_to_leaf = reverse(leaf_to_root)

    println("constuctZ")
    Z = ConstructZ(tree)
    (K, K) = size(model.weights)
    old_relevant_pairs = Array(Array{Int64, 2}, (K,K))
    println("construct relevant pairs")
    for k1 = 1:K
        for k2 = 1:K
            old_relevant_pairs[k1,k2] = compute_relevant_pairs(Z,k1,k2)
        end
    end 


    for node_index = root_to_leaf

        if rand() < .95
            continue
        end

        println("vardim sampling node_index: ", node_index)
        u = tree.nodes[node_index].state
        effective_lambda = model.lambda * (1 - model.gamma) * model.gamma^(tree.nodes[node_index].num_ancestors - 1)

        (K, K) = size(model.weights)
        
        W_index_pointers = weight_index_pointers(model.tree)
        start_index = W_index_pointers[node_index]
        end_index = start_index + u - 1
        new_k  = end_index + 1

        Z = ConstructZ(tree)
        # construct new Z col
        #println("construct Znew")
        leaves = GetLeaves(tree, node_index)
        Zk = zeros(Int64, N)
        Zk[leaves] = 1

        Znew = zeros(Int64, (N, K+1))
        Znew[:,1:end_index] = Z[:,1:end_index]
        Znew[:,end_index+1] = Zk
        Znew[:,end_index+2:end] = Z[:,end_index+1:end]        
        Znew = sparse(Znew)

        new_relevant_pairs = Array(Array{Int64, 2}, (K+1,K+1))
        # Leave space for new feature's weights
        #println("construct newmodel weights")
        new_model = copy(model)
        new_model.weights = zeros(Float64, (K+1,K+1))
        nonzero_element_indices = [x <= end_index ? x : x + 1 for x in 1:K]
        new_model.weights[nonzero_element_indices, nonzero_element_indices] = copy(model.weights)
        W = new_model.weights

        new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = copy(old_relevant_pairs)

        if model_spec.debug
            println("W[start,start]: ", new_model.weights[start_index,start_index])
        end

#        for k1 = nonzero_element_indices
#        new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = relevant_pairs

        #println("construct new relevant pairs")
        # update end_index+1 row and col of relevant pairs
        for k1 = 1:K+1
            k2 = new_k
            new_relevant_pairs[k1,k2] = compute_relevant_pairs(Znew,k1,k2)
            new_relevant_pairs[k2,k1] = compute_relevant_pairs(Znew,k2,k1)
        end 

        # All new weights initialized to zero, no need to update latent effects for the new feature
        # need to update them for the cases where we are removing features
        #println("construct component_latent_effects")
        component_latent_effects = Array(Array{Float64,2}, u+2)
        num_local_mutations = zeros(Int64, u+2)
        for u_ind = 1:u+2
            component_latent_effects[u_ind] = copy(latent_effects)
            if u_ind <= u
                num_local_mutations[u_ind] = u - 1
                removed_k = start_index + u_ind - 1
                for k = 1:K+1
                    for p = 1:size(new_relevant_pairs[removed_k, k])[2]
                        i = new_relevant_pairs[removed_k, k][1,p]
                        j = new_relevant_pairs[removed_k, k][2,p]
                        component_latent_effects[u_ind][i,j] -= W[removed_k, k]
                    end
                    if k != removed_k
                        for p = 1:size(new_relevant_pairs[k, removed_k])[2]
                            i = new_relevant_pairs[k, removed_k][1,p]
                            j = new_relevant_pairs[k, removed_k][2,p]
                            component_latent_effects[u_ind][i,j] -= W[k, removed_k]
                        end
                    end
                end
            elseif u_ind == u+1
                num_local_mutations[u_ind] = u
            else
                num_local_mutations[u_ind] = u + 1
            end
        end
   
        # Construct w_is_auxiliary for each component
        w_is_auxiliary = Array(Array{Bool,2},0)
        U = [ x <= u ? u-1 : x == u+1 ? u : u+1 for x = 1:u+2] # [u-1, u-1, u-1, ..., u-1, u, u+1]
        for u_ind = 1:u
            k_ind = start_index + u_ind - 1
            push(w_is_auxiliary,[x == k_ind  || y == k_ind ||
                                 x == new_k || y == new_k
                                 for x = 1:K+1, y = 1:K+1] )
        end

        push(w_is_auxiliary,[x === new_k || y == new_k for x = 1:K+1, y = 1:K+1] )
        push(w_is_auxiliary, zeros(Bool, (K+1,K+1)))

        # "constant" terms needed so that the sum of the different mixture components
        # is done with each scaled correctly
        const_terms = vardim_splits(new_model, Y, component_latent_effects,
                          observed_effects, U, W_index_pointers, node_index,
                          effective_lambda, w_is_auxiliary) 

        #println("slice sample")
        for iter = 1:num_W_sweeps
            for k1 = 1:K+1
                if k1 >= start_index && k1 <= end_index+1
                    k2_range = 1:K+1
                else
                    k2_range = start_index:end_index+1
                end
                for k2 = k2_range

                    #println("slice sample prep, start, end:", start_index, " ", end_index)

                    w_is_auxiliary_k1k2 = zeros(Bool, length(num_local_mutations))

                    aux1 = k1 - start_index + 1
                    aux2 = k2 - start_index + 1

                    if start_index <= end_index
                        if aux1 > 0 && aux1 <= u
                            w_is_auxiliary_k1k2[aux1] = true
                        end
                        if aux2 > 0 && aux2 <= u
                            w_is_auxiliary_k1k2[aux2] = true
                        end
                    end

                    #println("length(num_local_mutations): ", length(num_local_mutations))
                    if k1 == end_index+1 || k2 == end_index+1
                        w_is_auxiliary_k1k2 = ones(Bool, length(num_local_mutations))
                        w_is_auxiliary_k1k2[end] = false
                    end

                    for uuu = 1:length(w_is_auxiliary_k1k2)
                        assert(w_is_auxiliary_k1k2[uuu] == w_is_auxiliary[uuu][k1,k2])
                    end

                    w_cur = W[k1,k2]
                    w_old = w_cur

                    delta_const_terms = vardim_local_splits(new_model, Y, 
                        new_relevant_pairs[k1,k2], component_latent_effects,
                        observed_effects, num_local_mutations, W_index_pointers,
                        node_index, effective_lambda, zeros(length(num_local_mutations)),
                        w_old, w_old, w_is_auxiliary_k1k2) 

                    const_terms -= delta_const_terms

                    g = x -> vardim_local_sum(new_model, Y, new_relevant_pairs[k1,k2],
                                 component_latent_effects, observed_effects,
                                 num_local_mutations, W_index_pointers, node_index,
                                 effective_lambda, const_terms, w_old, x,
                                 w_is_auxiliary_k1k2) 
                    gx0 = g(w_cur)
                    for slice_iter = 1:num_W_slice_steps
                        (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                    end
 
                    W[k1,k2] = w_cur

                    if model_spec.debug
                        delta_const_terms = vardim_local_splits(new_model, Y, 
                            new_relevant_pairs[k1,k2], component_latent_effects,
                            observed_effects, num_local_mutations, W_index_pointers,
                            node_index, effective_lambda, zeros(length(num_local_mutations)),
                            w_old, w_cur, w_is_auxiliary_k1k2)
                        println("w_old, w_new, prob: ", w_old,", ", w_cur, ", ", g(w_cur))
                        println(delta_const_terms)
                        println(const_terms)
                    end

                    delta_const_terms = vardim_local_splits(new_model, Y, 
                        new_relevant_pairs[k1,k2], component_latent_effects,
                        observed_effects, num_local_mutations, W_index_pointers,
                        node_index, effective_lambda, zeros(length(num_local_mutations)),
                        w_old, w_cur, w_is_auxiliary_k1k2) 

                    const_terms += delta_const_terms

                    #println("adjust precomputations")

                    for p = 1:size(new_relevant_pairs[k1,k2])[2]
                        i = new_relevant_pairs[k1,k2][1,p]
                        j = new_relevant_pairs[k1,k2][2,p]

                        latent_effects[i,j] += w_cur - w_old


                        for u_ind = find(!w_is_auxiliary_k1k2)
                            component_latent_effects[u_ind][i,j] += w_cur - w_old
                        end

                    end
                    
     
                end
            end
        end

        logprobs = vardim_splits(new_model, Y, component_latent_effects, observed_effects,
                                 U, W_index_pointers, node_index, effective_lambda,
                                 w_is_auxiliary)

        sampled_component = randmult(exp_normalize(logprobs)) 

        new_u = U[sampled_component]
        new_K = K + new_u - u
        new_valid_indices = find(!w_is_auxiliary[sampled_component])
        new_W = reshape(copy(W[new_valid_indices]), (new_K, new_K)) 
        old_relevant_pairs = reshape(new_relevant_pairs[new_valid_indices], (new_K, new_K))

        oldW = copy(model.weights)

        model.weights = copy(new_W)
        tree.nodes[node_index].state = new_u

#        tree_prior = prior(model) 
#        tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#        println("new probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

        latent_effects = component_latent_effects[sampled_component]

#        println("new_u,weight_prob,new_weight_prob ", new_u, ",", sum(normal_logpdf(oldW, model.w_sigma )), ",", sum(normal_logpdf(new_W, model.w_sigma )))

        println("K,eff_lambda,Delta_K: ", new_K, ",",effective_lambda, ",",new_u - u)
        println(logprobs)

        if model_spec.debug

            Z = ConstructZ(model.tree)

            full_latent_effect = Z * model.weights * Z'

            println("maxdiff latent effect: ", max(abs(full_latent_effect-latent_effects)))

            full_probs = zeros(length(logprobs))
            for cmp = 1:length(logprobs)

                new_uu = U[cmp]
                new_K = K + new_uu - u
                new_valid_indices = find(!w_is_auxiliary[cmp])
                new_W = reshape(copy(W[new_valid_indices]), (new_K, new_K)) 

                new_model.weights = copy(new_W)
                new_model.tree.nodes[node_index].state = new_uu

                tree_prior = prior(new_model) 
                tree_LL = likelihood(new_model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
                full_probs[cmp] = tree_prior + tree_LL
                if cmp == sampled_component
                    str = "*tree probability["
                else
                    str = "tree probability["
                end
                println(str,cmp,"], prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
            end
        
            model.tree.nodes[node_index].state = new_u
            println("sampled_component: ", sampled_component)

            println(logprobs - max(logprobs))
            println(full_probs - max(full_probs))
#            tree_prior = prior(model) 
#            tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#            println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)


#            println("start: ", start_index)
#            println("new_k: ", new_k)
#            println("model[start,start]: ", model.weights[start_index,start_index])
#            println("W[start,start]: ", model.weights[start_index,start_index])
#            println("new_model[start,start]: ", new_model.weights[start_index,start_index])
#            println("newW[start,start]: ", new_W[start_index,start_index])
#            println("oldW[start,start]: ", oldW[start_index,start_index])
#            if new_u == u
#                println("Wdiff: ", norm(new_W-model.weights))
##                println(oldW)
##                println(model.weights)
##                println(new_model.weights)
#            end
            #println(!w_is_auxiliary[sampled_component])
        end
    end
end


function sample_psi(model::ModelState,
                    model_spec::ModelSpecification,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2},
                    latent_effects::Array{Float64,2},
                    observed_effects::Array{Float64,2})

    tree = model.tree
    N::Int = (length(tree.nodes)+1)/2

#    tree_prior = prior(model) 
#    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    #for prune_index = 1:2N-1
    for prune_index = 1:N # only prune leaves, it's much faster

        old_model = copy(model)
        old_tree = old_model.tree
        old_W = old_model.weights

        parent = tree.nodes[prune_index].parent
        if parent == Nil()
            continue
        end

        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil()
            continue
        end

        println("Sampling Prune Index: ", prune_index, " Num Leaves: ", length(GetLeaves(tree, grandparent.index)))
        if parent.children[1].index == prune_index
            original_sibling = parent.children[2]
        else
            original_sibling = parent.children[1]
        end

        prune_tree!(model, prune_index)
        gp = grandparent.index


        println("Num Leaves pruned: ", length(GetLeaves(tree, prune_index)), " Num leaves remaining: ", length(GetLeaves(tree, gp)) )


        if rand() < model_spec.global_move_probability
            subtree_indices = GetSubtreeIndicies(tree, prune_index)
            i = 1
            while contains(subtree_indices, i)
                i += 1
            end
            root = FindRoot(tree, i)
            path = GetLeafToRootOrdering(tree, root.index)
        else
            leaf = GetRandomLeaf(tree, original_sibling.index)
            path = GetPath(tree, leaf)
        end

        sibling_path_index = find(path .== original_sibling.index)

        (priors, pstates) = psi_infsites_logpdf(model, prune_index, path)
        (likelihoods, lstates) = psi_likelihood_logpdf(model, model_spec, prune_index, path, Y, X_r, X_p, X_c)

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)

        if any(isnan(probs))
            nan_ind = find(isnan(probs))[1]
            maxprior = max(priors)
            maxll = max(likelihoods)
            println("ind,logprob,prob: ", (nan_ind, logprobs[nan_ind], probs[nan_ind]))
            println("maxprior, maxlikelihood = ", (maxprior, maxll))
            println(likelihoods)
            println(lstates)
            assert(false)
        end

        state_index = randmult(probs)

        subtree_indices = GetSubtreeIndicies(tree, prune_index)
        i = 1
        while contains(subtree_indices, i)
            i += 1
        end
        root = FindRoot(tree, i)

        A = Int64[]
        graft_indices = Int64[]
        new_parent_features = Array(Array{Int64,1},0)
        new_graftpoint_features = Array(Array{Int64,1},0)
        for i = 1:length(pstates)
            if pstates[i][1] == original_sibling.index
                push(A, i)
            end 

            push(graft_indices, pstates[i][1])
            push(new_graftpoint_features, pstates[i][2])
            push(new_parent_features, pstates[i][3])
        end



        assert(pstates[state_index] == lstates[state_index])

        (graft_index, graftpoint_features, parent_features) = pstates[state_index]

        graft_tree!(model, prune_index, graft_index, parent_features, graftpoint_features)


        if model_spec.debug
            println("original_index,insert_index,parent,root: ", original_sibling.index, ",", pstates[state_index][1], ",", parent.index, ",", root.index)
            println("logprobs: ", logprobs)
            println("graft indices: ", graft_indices)
            println("new parent features: ", new_parent_features)
            println("new graftpoint features: ", new_graftpoint_features)
            println("start_probs, end_prob: ",logprobs[A]  ,", ", logprobs[state_index])
            println("start_prior, start_LL: ",priors[A], ", ", likelihoods[A])
            println("end_prior, end_LL: ", priors[state_index], ", ", likelihoods[state_index])

            tree_prior = prior(model) 
            tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
            println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
        end



    end
end
# Integrate exp(f) in a numerically stable way
function int_exp(f::Function, a::Float64, b::Float64)
    fx::Float64 = f(0.5(a+b))::Float64
    g = (t::Float64) -> exp(f(t)::Float64-fx)
    (integral, error) = adaptive_gauss_kronrod(g,a,b,1e-10,50)

    (integral,fx)
end

function sample_time(f::Function,
                     a::Float64,
                     b::Float64,
                     num_iter::Int64,
                     step_scale::Float64)
    x = 0.5(a+b)
    fx = f(x)
    w = (b-a)*step_scale
    for i = 1:num_iter
        (x,fx) = slice_sampler(x, f, w, 10000, a, b, fx)
    end
    return x
end

#end #end profile
