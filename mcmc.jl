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

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prior = prior(model) 
        tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
        println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
        if iter == 2
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

    vardim_time = time()
    sample_Z(model, model_spec, Y, X_r, X_p, X_c, latent_effects, observed_effects)
    vardim_time = time() - vardim_time

    tree_prior = prior(model) 
    tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)


    psi_time = time()
    sample_psi(model, model_spec, Y, X_r, X_p, X_c, latent_effects, observed_effects)
    psi_time = time() - psi_time



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
                g = x -> W_local_logpdf(model, Y, relevant_pairs, latent_effects, observed_effects, w_cur, x)

                w_old = w_cur
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
            old_relevant_pairs[k2,k1] = old_relevant_pairs[k1,k2]
        end
    end 


    for node_index = root_to_leaf
        println("vardim sampling node_index: ", node_index)
        u = tree.nodes[node_index].state
        effective_lambda = model.lambda * (1 - model.gamma) * model.gamma^(tree.nodes[node_index].num_ancestors - 1)

        (K, K) = size(model.weights)
        
        W_index_pointers = weight_index_pointers(model.tree)
        start_index = W_index_pointers[node_index]
        end_index = start_index + u - 1


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
        new_model.weights[nonzero_element_indices, nonzero_element_indices] = model.weights
        W = new_model.weights

        new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = old_relevant_pairs

#        for k1 = nonzero_element_indices
#        new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = relevant_pairs

        #println("construct new relevant pairs")
        new_k  = end_index + 1
        # update end_index+1 row and col of relevant pairs
        for k1 = 1:K+1
            k2 = new_k
            new_relevant_pairs[k1,k2] = compute_relevant_pairs(Znew,k1,k2)
            new_relevant_pairs[k2,k1] = new_relevant_pairs[k1,k2]
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
    
        #println("slice sample")
        for iter = 1:num_W_sweeps
            for k1 = start_index:end_index+1
                for k2 = 1:K+1

                    #println("slice sample prep, start, end:", start_index, " ", end_index)

                    w_is_auxiliary = zeros(Bool, length(num_local_mutations))

                    aux1 = k1 - start_index + 1
                    aux2 = k2 - start_index + 1

                    if start_index <= end_index
                        w_is_auxiliary[aux1] = true
                        if aux2 > 0 && aux2 <= u
                            w_is_auxiliary[aux2] = true
                        end
                    end

                    #println("length(num_local_mutations): ", length(num_local_mutations))
                    if k1 == end_index+1 || k2 == end_index+1
                        w_is_auxiliary = ones(Bool, length(num_local_mutations))
                        w_is_auxiliary[end] = false
                    end



                    w_cur = W[k1,k2]
                    #println("sample")

                    w_old = w_cur
                    g = x -> vardim_local_sum(new_model, Y, new_relevant_pairs[k1,k2], component_latent_effects, observed_effects,
                                              num_local_mutations, W_index_pointers, node_index, effective_lambda, w_old,
                                              x, w_is_auxiliary) 
                    gx0 = g(w_cur)
                    for slice_iter = 1:num_W_slice_steps
                        (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                    end
                  
                    W[k1,k2] = w_cur


                    #println("adjust precomputations")

                    for p = 1:size(new_relevant_pairs[k1,k2])[2]
                        i = new_relevant_pairs[k1,k2][1,p]
                        j = new_relevant_pairs[k1,k2][2,p]

                        latent_effects[i,j] += w_cur - w_old


                        for u_ind = find(!w_is_auxiliary)
                            component_latent_effects[u_ind][i,j] += w_cur - w_old
                        end

                    end
                    
     
                end
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

        logprobs = vardim_splits(new_model, Y, component_latent_effects, observed_effects, U, W_index_pointers, 
                                 node_index, effective_lambda, w_is_auxiliary) 

        sampled_component = randmult(exp_normalize(logprobs)) 

        new_u = U[sampled_component]
        new_K = K + new_u - u
        new_valid_indices = find(!w_is_auxiliary[sampled_component])
        new_W = reshape(W[new_valid_indices], (new_K, new_K)) 
        old_relevant_pairs = reshape(new_relevant_pairs[new_valid_indices], (new_K, new_K))

        oldW = copy(model.weights)

        model.weights = new_W
        tree.nodes[node_index].state = new_u

#        tree_prior = prior(model) 
#        tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#        println("new probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

        latent_effects = component_latent_effects[sampled_component]

#        println("new_u,weight_prob,new_weight_prob ", new_u, ",", sum(normal_logpdf(oldW, model.w_sigma )), ",", sum(normal_logpdf(new_W, model.w_sigma )))

        println("K,eff_lambda,Delta_K: ", new_K, ",",effective_lambda, ",",new_u - u)
        println(logprobs)

        if model_spec.debug
            tree_prior = prior(model) 
            tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
            println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
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


        parent = tree.nodes[prune_index].parent
        if parent == Nil()
            continue
        end

        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil()
            continue
        end

        println("Sampling Prune Index: ", prune_index, " Num Leaves: ", length(GetLeaves(tree, grandparent.index)))
        prune_tree!(model, prune_index)
        gp = grandparent.index


        println("Num Leaves pruned: ", length(GetLeaves(tree, prune_index)), " Num leaves remaining: ", length(GetLeaves(tree, gp)) )

        original_sibling = grandparent.children[find( grandparent.children .!= Nil())][1]

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
        (likelihoods, lstates) = psi_likelihood_logpdf(model, prune_index, path, Y, X_r, X_p, X_c)

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
        for i = 1:length(pstates)
            if pstates[i][1] == original_sibling.index
                push(A, i)
            end 

            push(graft_indices, pstates[i][1])
            push(new_parent_features, pstates[i][3])
        end


#        println("original_index,insert_index,root: ", original_sibling.index, ",", pstates[state_index][1], ",", root.index)
#        println("logprobs: ", logprobs)
#        println("graft indices: ", graft_indices)
#        println("new parent features: ", new_parent_features)
#        println("start_probs, end_prob: ",logprobs[A]  ,", ", logprobs[state_index])
#        println("start_prior, start_LL: ",priors[A], ", ", likelihoods[A])
#        println("end_prior, end_LL: ", priors[state_index], ", ", likelihoods[state_index])


        assert(pstates[state_index] == lstates[state_index])

        (graft_index, graftpoint_features, parent_features) = pstates[state_index]

        graft_tree!(model, prune_index, graft_index, parent_features, graftpoint_features)

#        tree_prior = prior(model) 
#        tree_LL = likelihood(model, model_spec, Y, X_r, X_p, X_c, linspace(1,N,N))
#        println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)
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
