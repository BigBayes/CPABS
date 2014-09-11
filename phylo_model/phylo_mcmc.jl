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
            ZZ, Etas, inds = model2array(model)
            p_dendrogram = dendrogram(ZZ,zeros(length(inds)), plot=false)

            tbl = Table(1,1)
            tbl[1,1] = p_dendrogram
            #Winston.display(tbl)
 
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
    U = sum(u)

    Tau = Array(Node, 2N-1)
    P = Array(Array{Node}, 2N-1)
    #Self inclusive
    ancestors = Array(Array{Node}, 2N-1)
    K = zeros(2N-1, 2N-1)

    Pki = [Array(Node,0) for i = 1:2N-1]
    An_i = [Array(Node,0) for i = 1:2N-1]
    An_tau = [Array(Node,0) for i = 1:2N-1]

    A_I = zeros(2N-1)
    A_tau = zeros(2N-1)

    An_ni = [Array(Node,0) for i = 1:2N-1]
    An_ntau = [Array(Node,0) for i = 1:2N-1]

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

    for i = indices
        tau = Tau[i]

        An = tau == Nil() ? [] : ancestors[tau.index]

        for j in An
            U_i[i] += u[j.index]
        end 
    end

    println("blah1") 
    for i = indices
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
                if !(l in ancestors[i]) && !(r in ancestors[i])
                    push!(An_ni[k], tree.nodes[i])
                end
                tau = Tau[i]
                An = tau == Nil() ? [] : ancestors[tau.index]
                if !(l in An) && !(r in An)
                    push!(An_ntau[k], tree.nodes[i])
                end
            end 
        end
    end

    println("blah2") 

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


            for a in [l,r]
                # Compute K(a,n) for n \in {n | a \in P(n)} 
                for n = Pki[a]
                    j = n.index
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

            for n = An_ni[i]
                j = n.index
                prod = 1.0

                for k = ancestors[j]
                    if k == l || k == r
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
                B_I[i] += prod 
            end

            for n = An_ntau[i]
                j = n.index
                prod = 1.0
                An_k = Tau[j] == Nil() ? [] : ancestors[Tau[j].index]
                for k = An_k
                    if k == l || k == r
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
                B_tau[i] += prod 
            end
            
            C_l = A_tau[l] - A_I[l]
            C_r = A_tau[r] - A_I[r]
            D = B_tau[i] - B_I[i]

            p_s = 1 - 2/(N_l+N_r+1)

            # Sample nutd_l
            f = x -> logsumexp( nu_tilde_splits(nu_r, x, nutd_r, gam, U, U_i[l], U_i[r], u, 
                                    K[l,:], K[r,:], C_l, C_r, D, P[l], P[r], p_s, node="l"))



            nutd_l = nutd_l == 1.0 ? rand(Uniform(0,1)) : nutd_l
            (nutd_u, f_nutd) = slice_sampler(nutd_l, f, 0.1, 10, 0.0, 1.0)

            f_vals = nu_tilde_splits(nu_r, nutd_u, nutd_r, gam, U, U_i[l], U_i[r], u, 
                                    K[l,:], K[r,:], C_l, C_r, D, P[l], P[r], p_s, node="l")

            f_ind = rand(Categorical(exp_normalize(f_vals)))
            nutd_l = f_ind == 1 ? 1.0 : nutd_u

            # Sample nutd_r
            f = x -> logsumexp( nu_tilde_splits(nu_r, nutd_l, x, gam, U, U_i[l], U_i[r], u, 
                                    K[l,:], K[r,:], C_l, C_r, D, P[l], P[r], p_s, node="r"))

            nutd_r = nutd_r == 1.0 ? rand(Uniform(0,1)) : nutd_r
            (nutd_u, f_nutd) = slice_sampler(nutd_r, f, 0.1, 10, 0.0, 1.0)

            f_vals = nu_tilde_splits(nu_r, nutd_l, nutd_u, gam, U, U_i[l], U_i[r], u, 
                                    K[l,:], K[r,:], C_l, C_r, D, P[l], P[r], p_s, node="r")

            f_ind = rand(Categorical(exp_normalize(f_vals)))
            nutd_r = f_ind == 1 ? 1.0 : nutd_u

            # Sample nu_r = 1-nu_l

            f = x -> nu_logpdf(x, nutd_l, nutd_r, gam, U, U_i[l], U_i[r], u,
                                K[l,:], K[r,:], C_l, C_r, D, P[l], P[r], cur.state, alpha, N_l, N_r) 
            (nu_r, f_nu) = slice_sampler(nu_r, f, 0.1, 10, 0.0, 1.0)

            cur.rho = nu_r
            left_child.rhot = nutd_l
            right_child.rhot = nutd_r
        end
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


        correct_priors, correct_likelihoods = prune_graft_logprobs(model, model_spec, data, prune_index)
        correct_logprobs = correct_priors + correct_likelihoods

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

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)



        priors = priors .- maximum(priors)
        likelihoods = likelihoods .- maximum(likelihoods)
        logprobs = logprobs .- maximum(logprobs)

        correct_priors = correct_priors .- maximum(correct_priors)
        correct_likelihoods = correct_likelihoods .- maximum(correct_likelihoods)
        correct_logprobs = correct_logprobs - maximum(correct_logprobs)
       
        println("states: $pstates") 
        println("priors (efficient): $priors")
        println("priors (correct): $correct_priors")
        println("likelihoods (efficient): $likelihoods")
        println("likelihoods (correct): $correct_likelihoods")
        println("logprobs (efficient): $logprobs")
        println("logprobs (correct): $correct_logprobs")

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
        eta_logpdf(model, model_spec, data)
    end

    function eta_grad(eta::Vector{Float64})
        for j = N+1:2N-1
            tree.nodes[j].state = eta[1 + (j-N-1)*S : (j-N)*S]
        end
        eta_log_gradient(model, model_spec, data)
    end

    opts = @options numsteps=8 stepsize=0.1 transformation=ReducedNaturalTransformation

    for i = 1:hmc_iterations
        eta = hmc_sampler(eta, eta_density, eta_grad, opts)
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
