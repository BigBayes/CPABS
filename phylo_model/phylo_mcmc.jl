require("phylo_model/phylo_model.jl")
require("tree.jl")
require("utils/probability_util.jl")
require("phylo_model/pdf.jl")
require("samplers/transformation.jl")

plot_utils_loaded = true
if Pkg.installed("Winston") != nothing
    require("utils/plot_utils.jl")
else
    println("Failed to load plot_utils.jl: disabling plotting")
    plot_utils_loaded = false
end

require("samplers/slicesampler.jl")
require("samplers/refractive_sampler.jl")
require("samplers/hmc.jl")

#@profile begin
function mcmc(data::DataState,
              lambda::Float64,
              gam::Float64,
              alpha::Float64,
              init_N::Int64,
              model_spec::ModelSpecification,
              iterations::Int64,
              burnin_iterations::Int64)

    # number of leaves is the number of split points plus one
    N = init_K+1

    (M,S) = size(data.reference_counts)
    eta = [rand(S) for i = 1:2N-1]
    tree = Tree(eta)

    InitializeBetaSplits(tree, () -> rand(Beta(1,1)))


    # initial root node must have nutd < 1.0
    root = FindRoot(tree,1)
    tree.nodes[root.index].rhot = 0.9


    Z = rand(1:N-1, M) + N
    # choose one mutation at random for each cluster to ensure they are nonempty
    perm = randperm(M)
    for k = 1:N-1
        Z[perm[k]] = k + N
    end


    model = ModelState(lambda,gam,alpha,tree,Z)


    trainLLs = Float64[]
    Ks = Int[]

    iters = Int[]

    models = Array(ModelState,0)
    debug = model_spec.debug 
    model_spec.debug = false

    local tbl

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prior = prior(model,model_spec)
        tree_LL = likelihood(model, model_spec, data)
        println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

        if iter == 2 && debug
            model_spec.debug = true
        end

        mcmc_sweep(model, model_spec, data)

        if model_spec.plot && plot_utils_loaded

            u = zeros(Int64, 2N-1)
            for i=1:length(Z)
                u[Z[i]] += 1
            end 

            ZZ, leaf_times, Etas, inds = model2array(model, return_leaf_times=true)
            p_dendrogram = dendrogram(ZZ,u[inds], plot=false, sorted_inds=inds, leaf_times=leaf_times)

            tbl = Table(1,1)
            tbl[1,1] = p_dendrogram
            Winston.display(tbl)
 
        end


        if mod(iter, 10) == 0 || iter == iterations
            push!(models,copy(model))
            push!(iters, iter)
        end
        push!(trainLLs, tree_LL)
        #push!(Ks, size(model.weights)[1])
    end

    if model_spec.plot && plot_utils_loaded
        (iters, Ks, trainLLs, models, tbl )
    else
        (iters, Ks, trainLLs, models )
    end
end


function mcmc_sweep(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

    tree = model.tree

    tree_prior = prior(model,model_spec) 
    tree_LL = likelihood(model, model_spec, data)
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    psi_time = time()
    sample_psi(model, model_spec, data)
    psi_time = time() - psi_time

    tree_prior = prior(model,model_spec) 
    tree_LL = likelihood(model, model_spec, data)
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    slice_iterations = 5
    nu_time = time()
    sample_nu_nutd(model, model_spec, slice_iterations)
    nu_time = time() - nu_time
    tree_prior = prior(model,model_spec) 
    tree_LL = likelihood(model, model_spec, data)
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    hmc_iterations = 5
    eta_time = time()
    sample_eta(model, model_spec, data, hmc_iterations)
    eta_time = time() - eta_time

    tree_prior = prior(model,model_spec) 
    tree_LL = likelihood(model, model_spec, data)
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    Z_iterations = 1
    Z_time = time()
    sample_assignments(model, model_spec, data, Z_iterations)
    Z_time = time() - Z_time
    tree_prior = prior(model,model_spec) 
    tree_LL = likelihood(model, model_spec, data)
    println("tree probability, prior, LL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL)

    println("MCMC Timings (psi, nu, eta, Z) = ", (psi_time, nu_time, eta_time, Z_time))
end


# nu = nu_r = 1-nu_l = parent.rho
# nutd_l = l.rhot 
function sample_nu_nutd(model::ModelState,
                        model_spec::ModelSpecification,
                        slice_iterations::Int)

    println("Sample nu, nu-tilde")
    tree = model.tree
    gam = model.gamma
    lambda = model.lambda
    alpha = model.alpha
    Z = model.Z

    N::Int = (length(tree.nodes)+1)/2

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    u = zeros(2N-1)
    U_i = zeros(2N-1)

    Tau = Array(Node, 2N-1)
    P = Array(Array{Node}, 2N-1)
    #Self inclusive
    ancestors = Array(Array{Node}, 2N-1)
    K = zeros(2N-1, 2N-1)

    Pki = [Array(Node,0) for i = 1:2N-1] # {j | i in P[j]}
    An_i = [Array(Node,0) for i = 1:2N-1] # {j | i in ancestors[j]}
    An_tau = [Array(Node,0) for i = 1:2N-1] #{j | i in ancestor[tau(j)]}

    A_I = zeros(2N-1)
    A_tau = zeros(2N-1)

    An_ni = [Array(Node,0) for i = 1:2N-1] # {j | i not in ancestors[j]} 
    An_ntau = [Array(Node,0) for i = 1:2N-1] # {j | i not in ancestors[tau(j)]}

    root_An_ni = Array(Node,0)
    root_An_ntau = Array(Node,0)

    B_I = zeros(2N-1)
    B_tau = zeros(2N-1)

    C = [Int64[] for x=1:maximum(Z)]
    for i = 1:length(Z)
        push!(C[Z[i]], i)
    end

    for i = indices
        cur = tree.nodes[i]
        P[i] = tau_path(cur)
        Tau[i] = tau(cur)
        ancestors[i] = GetAncestors(tree,i)
        u[i] = length(C[i])-1
        @assert u[i] >= 0 || i <= N
    end
    u[1:N] = 0
    U = sum(u)

    for i = indices
        tau = Tau[i]

        An = tau == Nil() ? [] : ancestors[tau.index]
        for j in An
            U_i[j.index] += u[i]
        end 
    end

    
    for i = indices
        if i > N
            for k = P[i]
                push!(Pki[k.index], tree.nodes[i])
            end
            for k = ancestors[i]
                push!(An_i[k.index], tree.nodes[i])
            end
            tau = Tau[i]
            An = tau == Nil() ? [] : ancestors[tau.index]
            for k = An
                push!(An_tau[k.index], tree.nodes[i])
            end

            for k = indices
                cur = tree.nodes[k]
                if k > N
                    left_child = cur.children[2]
                    right_child = cur.children[1]
                    l = left_child.index
                    r = right_child.index
                   
     
                    if !(left_child in ancestors[i]) && !(right_child in ancestors[i])
                        push!(An_ni[k], tree.nodes[i])
                    end
                    tau = Tau[i]
                    An = tau == Nil() ? [] : ancestors[tau.index]
                    if !(left_child in An) && !(right_child in An)
                        push!(An_ntau[k], tree.nodes[i])
                    end
                end 
            end
           
            # (root_An_ni is always empty)

            if Tau[i] == Nil()
                push!(root_An_ntau, tree.nodes[i])
            end
        end
    end

    for i = reverse(indices)
        if i > N
            cur = tree.nodes[i]
            left_child = tree.nodes[i].children[2]
            right_child = tree.nodes[i].children[1]
            l = left_child.index
            r = right_child.index

            eta = cur.state

            N_l = left_child.num_leaves
            N_r = right_child.num_leaves
            N_p = cur.num_leaves
            
            nu_p = 1.0
            if i != root.index
                parent = cur.parent
                self_direction = find(parent.children .== cur)[1]
                nu_p = self_direction == 1 ? parent.rho : 1-parent.rho 
            end

            rhot = cur.rhot
            nutd_l = left_child.rhot
            nutd_r = right_child.rhot
 
            nu_r = cur.rho
            nu_l = 1-cur.rho

            update_nu_nutd_constants!(i,[l,r], gam, ancestors, Tau, Pki, P, An_i, An_tau, An_ni[i], An_ntau[i], K, A_I, A_tau, B_I, B_tau)

            K_l = zeros(2N-1)
            K_r = zeros(2N-1)


            times = compute_times(model)



            C_l = A_tau[l] - A_I[l]
            C_r = A_tau[r] - A_I[r]
            D = B_tau[i] - B_I[i]

#            ZZ, leaf_times, Etas, inds = model2array(model, return_leaf_times=true)
#            dend = dendrogram(ZZ,u[inds], leaf_times=leaf_times, sorted_inds=inds, plot=false)
#
#            v = zeros(N-1)
#            for j = 1:2N-1
#                jcur = tree.nodes[j]
#                parent = jcur.parent
#                if j > N
#                    tau_node = tau(jcur)
#                    tau_t = tau_node == Nil() ? 1.0 : times[tau_node.index]
#                    v[j-N] = tau_t - times[j]
#                end
#
#            end
#
#            V_fast = (nu_r*nutd_r)^gam * C_r + ((1-nu_r)*nutd_l)^gam * C_l + D 
#            V_slow = sum(v)
#
#            println("V_fast: $V_fast")
#            println("V_slow: $V_slow")
#            println("v: $v")
            p_s = 1 - 2/(N_l+N_r+1)
            p_sl = 1 - 2/(N_l+1)
            p_sr = 1 - 2/(N_r+1)

            

            if l > N

                # Sample nutd_l
                f = x -> logsumexp( nu_tilde_splits(nu_r, x, nutd_r, gam, U, U_i[l], U_i[r], u, 
                                        K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r], p_sl, node="l"))



                nutd_l = nutd_l == 1.0 ? rand(Uniform(0,1)) : nutd_l
                (nutd_u, f_nutd) = slice_sampler(nutd_l, f, 0.1, 10, 0.0, 1.0)

                f_vals = nu_tilde_splits(nu_r, nutd_u, nutd_r, gam, U, U_i[l], U_i[r], u, 
                                        K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r], p_sl, node="l")

                f_ind = rand(Categorical(exp_normalize(f_vals)))
                nutd_l = f_ind == 1 ? 1.0 : nutd_u
            end

            if r > N
                # Sample nutd_r
                f = x -> logsumexp( nu_tilde_splits(nu_r, nutd_l, x, gam, U, U_i[l], U_i[r], u, 
                                        K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r], p_sr, node="r"))

                nutd_r = nutd_r == 1.0 ? rand(Uniform(0,1)) : nutd_r
                (nutd_u, f_nutd) = slice_sampler(nutd_r, f, 0.1, 10, 0.0, 1.0)

                f_vals = nu_tilde_splits(nu_r, nutd_l, nutd_u, gam, U, U_i[l], U_i[r], u, 
                                        K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r], p_sr, node="r")

                f_ind = rand(Categorical(exp_normalize(f_vals)))
                nutd_r = f_ind == 1 ? 1.0 : nutd_u
            end

            # Sample nu_r = 1-nu_l

            f = x -> nu_logpdf(x, nutd_l, nutd_r, gam, U, U_i[l], U_i[r], u,
                                K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r], cur.state, alpha, N_l, N_r)

#            f1 = x -> p_z_given_nu_nutd(x, nutd_l, nutd_r, gam, U, U_i[l], U_i[r], u, 
#                                        K[l,:], K[r,:], C_l, C_r, D, Pki[l], Pki[r]) + p_nu_Nl_Nr(x, N_l, N_r)
#
#
#            left_child.rhot = nutd_l
#            right_child.rhot = nutd_r
#            f2 = x -> (tree.nodes[i].rho = x; full_pdf(model, model_spec, data))
#
#                    x_range = [0.01:0.01:0.99]
#                    f1_vals = [f(x) for x in x_range]
#                    f2_vals = [f2(x) for x in x_range]
#
#                    f1_vals -= maximum(f1_vals)
#                    f2_vals -= maximum(f2_vals)
#
#                    c1 = Curve(x_range, f1_vals, color="blue")
#                    c2 = Curve(x_range, f2_vals, color="red")
#
#                    p = FramedPlot()
#                    add(p, c2) 
#                    add(p, c1) 
#                    tbl = Table(1,2)
#                    tbl[1,1] = p
#                    tbl[1,2] = dend
#                    display(tbl)
#                    @assert false
 
            (nu_r, f_nu) = slice_sampler(nu_r, f, 0.1, 10, 0.0, 1.0)

            cur.rho = nu_r

            if i == root.index
                update_nu_nutd_constants!(i,[i], gam, ancestors, Tau, Pki, P, An_i, An_tau, root_An_ni, root_An_ntau, K, A_I, A_tau, B_I, B_tau)

                C_i = A_tau[i] - A_I[i]
                D = B_tau[i] - B_I[i]

                nutd_p = cur.rhot 
                nutd_p = nutd_p == 1.0 ? rand(Uniform(0,1)) : nutd_p

                f = x -> logsumexp( root_nu_tilde_splits(x, gam, U, U_i[i], u, 
                                        K[i,:], C_i, D, Pki[i], p_s))
                ##########
#                    println("i: $i")
#                    f1 = x -> root_nu_tilde_splits(x, gam, U, U_i[i], u, 
#                                            K[i,:], C_i, D, Pki[i], p_s)[2]
#                    f2 = x -> (tree.nodes[i].rhot = x; prior(model, model_spec, debug=false))
#
#                    x_range = [0.01:0.01:0.99]
#                    f1_vals = [f1(x) for x in x_range]
#                    f2_vals = [f2(x) for x in x_range]
#
#                    f1_vals -= maximum(f1_vals)
#                    f2_vals -= maximum(f2_vals)
#
#                    c1 = Curve(x_range, f1_vals, color="blue")
#                    c2 = Curve(x_range, f2_vals, color="red")
#
#                    p = FramedPlot()
#                    add(p, c2) 
#                    add(p, c1) 
#                    tbl = Table(1,2)
#                    tbl[1,1] = p
#                    tbl[1,2] = dend
#                    display(tbl)
#                    @assert false
                ##########
            
                (nutd_u, f_nutd) = slice_sampler(nutd_p, f, 0.1, 10, 0.0, 1.0)
                f_vals = root_nu_tilde_splits(nutd_u, gam, U, U_i[i], u,
                                        K[i,:], C_i, D, Pki[i], p_s)
                f_ind = rand(Categorical(exp_normalize(f_vals)))
                nutd_p = f_ind == 1 ? 1.0 : nutd_u
                @assert nutd_p != 1.0
                cur.rhot = nutd_p
            end


        end
    end

end

function update_nu_nutd_constants!(i, update_set, gam, ancestors, Tau, Pki, P, An_i, An_tau, An_ni, An_ntau, K, A_I, A_tau, B_I, B_tau)

    for a in update_set
        # Compute K(a,n) for n \in {n | a \in P(n)} 
        for n = Pki[a]
            j = n.index
            K[a,j] = 1.0
            for k = P[j]
                if k.index == a
                    continue
                end
                if k.parent == Nil()
                    K[a,j] *= (k.rhot)^gam
                else
                    k_direction = find(k.parent.children .== k)[1]
                    nu_k = k_direction == 1 ? k.parent.rho : 1-k.parent.rho
                    K[a,j] *= (k.rhot*nu_k)^gam
                end

            end
        end

        # Compute A(a,I)
        A_I[a] = 0.0
        for n = An_i[a]
            j = n.index
            prod = 1.0
            for k = ancestors[j]
                if k.index == a
                    continue
                end
                if k.parent == Nil()
                    prod *= k.rhot^gam
                else
                    k_direction = find(k.parent.children .== k)[1]
                    nu_k = k_direction == 1 ? k.parent.rho : 1-k.parent.rho
                    prod *= (k.rhot*nu_k)^gam
                end
            end
            A_I[a] += prod
        end
        # Compute A(a,tau)
        A_tau[a] = 0.0
        for n = An_tau[a]
            j = n.index
            prod = 1.0
            An_k = Tau[j] == Nil() ? [] : ancestors[Tau[j].index]
            for k = An_k
                if k.index == a
                    continue
                end
                if k.parent == Nil()
                    prod *= k.rhot^gam
                else
                    k_direction = find(k.parent.children .== k)[1]
                    nu_k = k_direction == 1 ? k.parent.rho : 1-k.parent.rho
                    prod *= (k.rhot*nu_k)^gam
                end
            end
            A_tau[a] += prod
        end

    end

    
    B_I[i] = 0.0
    for n = An_ni
        j = n.index
        prod = 1.0

        for k = ancestors[j]
            @assert !(k in update_set)
            if k.parent == Nil()
                prod *= k.rhot^gam
            else
                k_direction = find(k.parent.children .== k)[1]
                nu_k = k_direction == 1 ? k.parent.rho : 1-k.parent.rho
                prod *= (k.rhot*nu_k)^gam
            end
        end
        B_I[i] += prod 
    end

    B_tau[i] = 0.0
    for n = An_ntau
        j = n.index
        prod = 1.0
        An_k = Tau[j] == Nil() ? [] : ancestors[Tau[j].index]
        for k = An_k
            @assert !(k in update_set)
            if k.parent == Nil()
                prod *= k.rhot^gam
            else
                k_direction = find(k.parent.children .== k)[1]
                nu_k = k_direction == 1 ? k.parent.rho : 1-k.parent.rho
                prod *= (k.rhot*nu_k)^gam
            end
        end
        B_tau[i] += prod 
    end
end

function sample_psi(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

    tree = model.tree
    N::Int = (length(tree.nodes)+1)/2


    println("psi: ")
    for parent_prune_index = N+1:2N-1
         
        parent = tree.nodes[parent_prune_index]
        prune_index = parent.children[1].index

        if model_spec.debug
            old_prior = prior(model, model_spec)
            old_LL = likelihood(model, model_spec, data)
        end

        if mod(parent_prune_index-N+1,ceil(N/10)) == 0
            percent = ceil((parent_prune_index-N+1)/ceil(N/10))*10
            println(" ",percent , "% ")
        end

        old_model = copy(model)
        old_tree = old_model.tree


        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil()
            continue
        end

        if parent.children[1].index == prune_index
            original_sibling = parent.children[2]
        else
            original_sibling = parent.children[1]
        end


        #correct_priors, correct_likelihoods = prune_graft_logprobs(model, model_spec, data, prune_index)

        PruneIndexFromTree!(model.tree, prune_index)

        gp = grandparent.index

        subtree_indices = GetSubtreeIndicies(tree, prune_index)
        i = 1
        while in(i,subtree_indices)
            i += 1
        end
        root = FindRoot(tree, i)
        path = GetLeafToRootOrdering(tree, root.index)
      

        (priors, pstates) = psi_infsites_logpdf(model, data, prune_index, path)
        (likelihoods, lstates) = psi_observation_logpdf(model, model_spec, data, prune_index, path)

#        priors = priors .- maximum(priors)
#
#        correct_priors = correct_priors .- maximum(correct_priors)
#
#        println("prune_index: $prune_index") 
#        println("parent_prune_index: $parent_prune_index") 
#        println("states: $lstates") 
#        println("priors (efficient): $priors")
#        println("priors (correct): $correct_priors")

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)

        if any(isnan(probs))
            nan_ind = find(isnan(probs))[1]
            maxprior = maximum(priors)
            maxll = maximum(likelihoods)
            println("ind,logprob,prob: ", (nan_ind, logprobs[nan_ind], probs[nan_ind]))
            println("maxprior, maxlikelihood = ", (maxprior, maxll))
            println(priors)
            println(likelihoods)
            println(lstates)
            assert(false)
        end

        state_index = randmult(probs)

        subtree_indices = GetSubtreeIndicies(tree, prune_index)
        i = 1
        while in(i, subtree_indices)
            i += 1
        end
        root = FindRoot(tree, i)



        @assert pstates[state_index] == lstates[state_index]

        graft_index = pstates[state_index]

        InsertIndexIntoTree!(model.tree, prune_index, graft_index) 

        println("graft_index: $graft_index")
        println("original_sibling_index: $(original_sibling.index)")

        if model_spec.debug 
            println("Sampling Prune Index: ", prune_index, " Num Leaves: ", length(GetLeaves(tree, grandparent.index)))
            println("Num Leaves pruned: ", length(GetLeaves(tree, prune_index)), " Num leaves remaining: ", length(GetLeaves(tree, gp)) )
            println("original_index,insert_index,parent,root: ", original_sibling.index, ",", pstates[state_index][1], ",", parent.index, ",", root.index)
#            println("logprobs: ", logprobs)
#            println("graft indices: ", graft_indices)

            println("graftpoint_features: $graftpoint_features")
            println("parent_features: $parent_features")

            subtree_indices = GetLeafToRootOrdering(tree, prune_index)
            grafttree_indices = GetLeafToRootOrdering(tree, graft_index)
            ancestor_indices = GetPath(tree, tree.nodes[prune_index].parent.index)

            count_features = x -> sum([tree.nodes[i].state for i in x])
            get_features = x -> [tree.nodes[i].state for i in x]

            subtree_num = count_features(subtree_indices)
            grafttree_num = count_features(grafttree_indices)
            ancestor_num = count_features(ancestor_indices)
           
            println("subtree num_features: $subtree_num") 
            println("grafttree num_features: $grafttree_num") 
            println("ancestor num_features: $ancestor_num") 
            println("original sibling under graftpoint?: $(original_sibling.index in grafttree_indices)")

            println("ancestors: $ancestor_indices")
            println("ancestor_features: $(get_features(ancestor_indices))")

            println("local_LL: $(likelihoods[state_index])")
            println("local_prior: $(priors[state_index])")
            println("old_local_LL: $(likelihoods[A])")
            println("old_local_prior: $(priors[A])")

            println("prob: $(probs[state_index])")
            println("old_prob: $(probs[A])")


            tree_prior = prior(model, model_spec) 
            tree_LL = likelihood(model, model_spec, data)
            println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

            full_diff = tree_LL + tree_prior - old_LL - old_prior
            local_diff = likelihoods[state_index] + priors[state_index] -
                         likelihoods[A[1]] - priors[A[1]]

            all_local_diff =  likelihoods[state_index] + priors[state_index] .- likelihoods[A] .- priors[A]

            min_local_diff, ai = findmin(abs(full_diff .- all_local_diff))

            #if (tree_LL-old_LL) + (tree_prior-old_prior) < -8
            if min_local_diff > 0.1 
                println("full_diff: $full_diff")
                println("local_diff: $local_diff")
                println("min_local_diff: $min_local_diff")

                println("prior_diff: $((tree_prior-old_prior))")
                println("local_prior_diff: $(priors[state_index] - priors[A[ai]])")

                println("prior_err: $(abs((tree_prior-old_prior) - (priors[state_index] - priors[A[1]])))")

                println("LL_diff: $(tree_LL-old_LL)")
                println("local_LL_diff: $(likelihoods[state_index] - likelihoods[A[ai]])")
                assert(false)
            end
        end
    end
end



function sample_eta(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState,
                    hmc_iterations::Int64)

    tree  = model.tree
    N::Int = (length(tree.nodes) + 1) / 2
    M, S = size(data.reference_counts)
    eta = zeros(S*(N-1))

    num_iterations = 10

    for j = N+1:2N-1
        eta[1 + (j-N-1)*S : (j-N)*S] = tree.nodes[j].state
    end

    function eta_density(eta::Vector{Float64})
        for j = N+1:2N-1
            tree.nodes[j].state = eta[1 + (j-N-1)*S : (j-N)*S]
        end
        e_pdf = eta_logpdf(model, model_spec, data)
        return e_pdf
    end

    function eta_grad(eta::Vector{Float64})
        for j = N+1:2N-1
            tree.nodes[j].state = eta[1 + (j-N-1)*S : (j-N)*S]
        end
        grad = eta_log_gradient(model, model_spec, data)
        return grad
    end

    hmcopts = @options numsteps=8 stepsize=0.1 transformation=ReducedNaturalTransformation
    refopts = @options m=2 w=0.5 refractive_index_ratio=1.3 transformation=ReducedNaturalTransformation

    for i = 1:hmc_iterations
        #eta = hmc_sampler(eta, eta_density, eta_grad, hmcopts)
        eta = refractive_sampler(eta, eta_density, eta_grad, refopts)
    end

    
 
    for j = N+1:2N-1
        tree.nodes[j].state = eta[1 + (j-N-1)*S : (j-N)*S]
    end
end

function sample_assignments(model::ModelState,
                            model_spec::ModelSpecification,
                            data::DataState,
                            num_iterations::Int64)

    tree = model.tree
    N::Int = (length(tree.nodes) + 1) / 2
    M,S = size(data.reference_counts)

    t = compute_times(model)
    Tau = compute_taus(model)
    phi = compute_phis(model)

    Z = model.Z
    U = zeros(Int64, N-1)
    for i=1:length(Z)
        U[Z[i]-N] += 1
    end 
    U = U .- 1

    for iter = 1:num_iterations    
        for i = 1:M
            cur_z = Z[i]-N
            z_probs = z_logpdf(model, model_spec, data, i, U, t, Tau, phi)

#            full_logpdf = k -> (model.Z[i] = k+N; full_pdf(model, model_spec, data))
#            full_z_probs = [full_logpdf(k) for k = 1:N-1]
#            z_probs -= maximum(z_probs)
#            full_z_probs -= maximum(full_z_probs)
#            println("full_z_probs: $full_z_probs")
#            println("z_probs: $z_probs")
#            @assert false


            new_z = rand(Categorical(exp_normalize(z_probs)))
            U[cur_z] -= 1
            U[new_z] += 1
            Z[i] = new_z+N
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
