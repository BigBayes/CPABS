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
        if u > .95
            U[i] = 2
        elseif u > .8
            U[i] = 1
        end
    end

    W = randn(sum(U),sum(U))
    tree = Tree(U)

    model = ModelState(lambda,gamma,1.0,1.0,tree,W,[0.0],[0.0],[0.0],[0.0],[0.0],0.0)

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prob = full_pdf(model, model_spec, Y, X_r, X_p, X_c)
        println("Tree probability: ", tree_prob)
        mcmc_sweep(model, model_spec, Y, X_r, X_p, X_c)
    end
end


function mcmc_sweep(model::ModelState,
                    model_spec::ModelSpecification,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2})

    tree = model.tree
    N::Int = (length(tree.nodes)+1)/2

    psi_time = time()

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


        leaf = GetRandomLeaf(tree, gp)
        path = GetPath(tree, leaf)

        (priors, pstates) = psi_infsites_logpdf(model, prune_index, path[1:end-1])
        (likelihoods, lstates) = psi_likelihood_logpdf(model, prune_index, path[1:end-1], Y, X_r, X_p, X_c)

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)

        if any(isnan(probs))
            nan_ind = find(isnan(probs))[1]
            maxprior = max(priors)
            maxll = max(likelihoods)
            println("ind,logprob,prob: ", (nan_ind, logprobs[nan_ind], probs[nan_ind]))
            println("maxprior, maxlikelihood = ", (maxprior, maxll))
            assert(false)
        end

        state_index = randmult(probs)

        (graft_index, graftpoint_features, parent_features) = pstates[state_index]

        graft_tree!(model, prune_index, graft_index, parent_features, graftpoint_features)
    end

    psi_time = time() - psi_time

    W_time = time()
    (K,K) = size(model.weights)

    num_W_sweeps = 1
    num_W_slice_steps = 1
    Z = ConstructZ(tree)
    W = model.weights

    latent_effects = zeros(Float64, (N,N))
    observed_effects = zeros(Float64, (N,N))
   
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c
 
    for i = 1:N
        for j = 1:N
            latent_effects[i,j] = (Z[i,:] * W * Z[j,:]')[1]
            observed_effects[i,j] = compute_observed_effects(model, model_spec, i, j, X_r, X_p, X_c)
        end
    end

    for iter = 1:num_W_sweeps
   
        #for k1 = 1:3
        for k1 = 1:K
            println("W sampling k1: ", k1)
            #for k2 = 1:3
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
            end

        end 
    end

    W_time = time() - W_time

    # Sample new features from root to leaf
    root = FindRoot(tree, 1) 
    leaf_to_root = GetLeafToRootOrdering(tree, root.index)
    root_to_leaf = reverse(leaf_to_root)


    vardim_time = time()
    for node_index = root_to_leaf
        println("vardim sampling node_index: ", node_index)
        u = tree.nodes[node_index].state
        effective_lambda = model.lambda * (1 - model.gamma) * model.gamma^(tree.nodes[node_index].num_ancestors - 1)

        (K, K) = size(model.weights)
        new_model = copy(model)
        new_model.weights = zeros(Float64, (K+1,K+1))
        
        W_index_pointers = weight_index_pointers(model.tree)
        start_index = W_index_pointers[node_index]
        end_index = start_index + u - 1


        Z = ConstructZ(tree)
        # construct new Z col
        println("construct Znew")
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
        println("construct newmodel weights")
        nonzero_element_indices = [x <= end_index ? x : x + 1 for x in 1:K]
        new_model.weights[nonzero_element_indices, nonzero_element_indices] = model.weights
        W = new_model.weights

#        for k1 = nonzero_element_indices
#        new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = relevant_pairs

        new_k  = end_index + 1
        # update end_index+1 row and col of relevant pairs
        for k1 = 1:K+1
            for k2 = 1:K+1
                new_relevant_pairs[k1,k2] = compute_relevant_pairs(Znew,k1,k2)
                new_relevant_pairs[k2,k1] = compute_relevant_pairs(Znew,k2,k1)
            end
        end 

        # All new weights initialized to zero, no need to update latent effects for the new feature
        # need to update them for the cases where we are removing features
        println("construct component_latent_effects")
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
    
        println("construct w norms") 
        w_k_norm = zeros(Float64, u+2)
        for u_ind = 1:u
            removed_k = start_index + u_ind - 1
            new_k = end_index+1

            w_k_norm[u_ind] = sum(W[removed_k,:].^2) + sum(W[:,removed_k].^2) -
                              W[removed_k, new_k]^2 - W[new_k, removed_k]^2 - W[removed_k, removed_k]^2
        end
        w_kp1_norm = 0.0
 
        println("slice sample") 
        for k1 = start_index:end_index+1
            for k2 = 1:K+1

                println("slice sample prep, start, end:", start_index, " ", end_index)

                w_is_auxiliary = zeros(Bool, length(num_local_mutations))

                aux1 = k1 - start_index + 1
                aux2 = k2 - start_index + 1

                if start_index <= end_index
                    w_is_auxiliary[aux1] = true
                    if aux2 > 0 && aux2 <= u
                        w_is_auxiliary[aux2] = true
                    end
                end

                println("length(num_local_mutations): ", length(num_local_mutations))
                if k1 == end_index+1 || k2 == end_index+1
                    w_is_auxiliary = ones(Bool, length(num_local_mutations))
                    w_is_auxiliary[end] = false
                end



                w_cur = W[k1,k2]
                println("sample")

                w_old = w_cur
                g = x -> vardim_local_sum(model, Y, new_relevant_pairs[k1,k2], component_latent_effects, observed_effects,
                                          num_local_mutations, W_index_pointers, node_index, effective_lambda, w_old,
                                          x, w_is_auxiliary) 
                gx0 = g(w_cur)
                for slice_iter = 1:num_W_slice_steps
                    (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                end
              
                W[k1,k2] = w_cur


                println("adjust precomputations")

                for p = 1:size(new_relevant_pairs[k1,k2])[2]
                    i = new_relevant_pairs[k1,k2][1,p]
                    j = new_relevant_pairs[k1,k2][2,p]

                    latent_effects[i,j] += w_cur - w_old


                    for u_ind = find(!w_is_auxiliary)
                        component_latent_effects[u_ind][i,j] += w_cur - w_old
                    end

                end
                
                for u_ind = find(w_is_auxiliary)
                    w_k_norm[u_ind] += w_cur^2 - w_old^2 
                end

                w_kp1_norm += w_cur^2 - w_old^2
 
            end
        end



    end

    vardim_time = time() - vardim_time

    println("MCMC Timings (psi, W, vardim) = ", (psi_time, W_time, vardim_time))
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
