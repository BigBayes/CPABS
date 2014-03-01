require("model.jl")
require("tree.jl")
require("probability_util.jl")
require("pdf.jl")

require("slicesampler.jl")
require("refractive_sampler.jl")
require("hmc.jl")

#@profile begin
function mcmc(data::DataState,
              lambda::Float64,
              gam::Float64,
              w_sigma::Float64,
              model_spec::ModelSpecification,
              iterations::Int64,
              burnin_iterations::Int64)

    YY = data.Ytrain
    X_r = data.X_r
    X_p = data.X_p
    X_c = data.X_c
 
    (N,N) = size(YY[1])

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
    W = randn(0,0)
    Waug = AugmentedMatrix(2N-1)
    tree = Tree(U)
    InitializeBetaSplits(tree, () -> rand(Beta(1,1)))
    model = ModelState(lambda,gam,w_sigma,1.0,w_sigma,tree,W,Waug,[0.0],[0.0],[0.0],[0.0],[0.0],0.0)

    for i = 1:2N-2 # the root should have no features
        num_ancestors = tree.nodes[i].num_ancestors
        #hard code lambda for initialization so we always start with some features
        effective_lambda = 4 * (0.5)^(num_ancestors-1)
        U[i] = effective_lambda > 10.0^-5 ? randpois(effective_lambda) : 0

        if U[i] > 0
            U[i] = 1
        end
        tree.nodes[i].state = U[i] 
    end

    sU = sum(U)
    if model_spec.diagonal_W
        w = randn(sU)
        W = diagm(w)
    else
        W = randn(sU,sU)
    end
    model.weights = W
    model.augmented_weights = AugmentedMatrix(2N-1, W, U)


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
#    parent_features = [start_index + st1:end_index] 
#    sibling_features = [start_index:start_index + st1 - 1]
#    graft_tree!(model,prune_index, sib_ind, parent_features, sibling_features)
#
#    assert(tree_copy == tree)
#    assert(norm(W) == norm(model.weights))
#    assert(W == model.weights)

    ###

    # if not diagonal, initialize using a diagonal model
    diagonal = model_spec.diagonal_W
#    model_spec.diagonal_W = true

    avg_test_likelihoods = [zeros(Float64, 0) for x in 1:N, y in 1:N]
    logit_args = [zeros(Float64, 0) for x = 1:N, y = 1:N]

    testLLs = Float64[]
    trainLLs = Float64[]
    Ks = Int[]

    train_errors = Float64[]
    test_errors = Float64[]
    avg_test_LLs = Float64[]
    aucs = Float64[]

    iters = Int[]

    models = Array(ModelState,0)
    debug = model_spec.debug 
    model_spec.debug = false

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prior = prior(model,model_spec) 
        tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
        println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

        if iter == 2 && debug
            model_spec.debug = true
        end
        if iter == 1
            model_spec.diagonal_W = diagonal
        end

        mcmc_sweep(model, model_spec, data)


        if iter > burnin_iterations
            Z = ConstructZ(model.tree) 
            testI, testJ = findn(data.Ytest[1] .>= 0)
            total_LL = 0.0

            for i = 1:N
                for j = 1:N
                    testLL,logit_arg = test_likelihood_ij(model,model_spec,data,Z,i,j)
                    push!(logit_args[i,j], logit_arg)

                    if data.Ytest[1][i,j] >= 0
                        push!(avg_test_likelihoods[i,j], testLL)
                        total_LL += logsumexp(avg_test_likelihoods[i,j]) - log(length(avg_test_likelihoods[i,j]))
                    end
                end
            end

            println("Average Test LL: ", total_LL / (length(testI)) )
            println("Total Test LL: ", total_LL  )

        end

        if mod(iter, 10) == 0 || iter == iterations
            push!(models,copy(model))
            push!(iters, iter)
            if iter > burnin_iterations
                (train_error, test_error, auc) = error_and_auc(logit_args, data)
                println("Train, Test, AUC: ", (train_error, test_error, auc))

                push!(train_errors, train_error)
                push!(test_errors, test_error)
                push!(avg_test_LLs, total_LL)
                push!(aucs, auc)
            end
        end
        push!(trainLLs, tree_LL)
        push!(testLLs, test_LL)
        push!(Ks, size(model.weights)[1])
    end

    (iters, train_errors, test_errors, avg_test_LLs, aucs, Ks, trainLLs, testLLs, models )
end


function mcmc_sweep(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

    YY = data.Ytrain
    (N,N) = size(YY[1])
    tree = model.tree
    Z = ConstructZ(tree)

    (latent_effects, observed_effects) = 
        construct_effects(model, model_spec, data, Z)

    sample_bias(model, model_spec, data, latent_effects, observed_effects)
    tree_prior = prior(model,model_spec) 
    tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
    println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

    psi_time = time()
    sample_psi(model, model_spec, data, latent_effects, observed_effects)
    psi_time = time() - psi_time

    tree_prior = prior(model,model_spec) 
    tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
    println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

    slice_iterations = 1
    rho_time = time()
    sample_rho_rhot(model, model_spec, slice_iterations)
    rho_time = time() - rho_time
    tree_prior = prior(model,model_spec) 
    tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
    println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

    W_time = time()
    Z = ConstructZ(tree)
    (latent_effects, observed_effects) = 
        construct_effects(model, model_spec, data, Z)

    #sample_W(model, model_spec, data, Z, latent_effects, observed_effects)
    sample_W_full(model, model_spec, data, Z, observed_effects)

    W_time = time() - W_time

    tree_prior = prior(model,model_spec) 
    tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
    println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)

    vardim_time = time()
    if rand() > 0.0
        sample_Z(model, model_spec, data, latent_effects, observed_effects,
                 tree_prior + tree_LL)
    end
    vardim_time = time() - vardim_time

    (K,K) = size(model.weights)
    println("MCMC Timings (psi, rho, W, vardim) = ", (psi_time, rho_time, W_time, vardim_time))
    println("K: ", K)
end

function sample_bias(model::ModelState,
                     model_spec::ModelSpecification,
                     data::DataState,
                     latent_effects::Array{Float64,2},
                     observed_effects::Array{Float64,2})

    println("Sample bias")
    num_slice_steps = 10

    c = model.c
    c_old = c
    g = x -> bias_logpdf(model,model_spec,data,latent_effects,observed_effects,
                         c_old, x)
    gx0 = g(c)

    for iter = 1:num_slice_steps
        (c, gx0) = slice_sampler(c, g, 1.0, 10, -Inf, Inf, gx0)
    end

    observed_effects += c - model.c
    model.c = c
end

function sample_W(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState,
                  Z,#sparse or full 2d array
                  latent_effects::Array{Float64,2},
                  observed_effects::Array{Float64,2})

    YY = data.Ytrain
    (K,K) = size(model.weights)
    (N,N) = size(YY[1])
    W = model.weights

    num_W_sweeps = 3
    num_W_slice_steps = 3

   
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_c

    Ynn = zeros(size(YY[1]))
    Ynn[find(YY[1] .>= 0)] = 1
    Ynn = sparse(Ynn)

    for iter = 1:num_W_sweeps
        println("W sampling sweep ", iter,"/",num_W_sweeps)



        for k1 = 1:K
            k2_range = model_spec.diagonal_W ? k1 :
                       model_spec.symmetric_W ? (1:k1) :
                       (1:K)
            for k2 = k2_range
                relevant_pairs = compute_relevant_pairs(model_spec,Ynn,Z,k1,k2)
                w_cur = W[k1,k2]
                w_old = w_cur
                g = x -> W_local_logpdf(model, data, relevant_pairs, latent_effects, observed_effects, w_old, x)


                gx0 = g(w_cur)
                gx00 = gx0
                for slice_iter = 1:num_W_slice_steps
                    (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                end
    
#                if model_spec.debug
#                    println(gx0-gx00)
#                end
 
                W[k1,k2] = w_cur

                
                for p = 1:size(relevant_pairs)[2]
                    i = relevant_pairs[1,p]
                    j = relevant_pairs[2,p]

                    latent_effects[i,j] += w_cur - w_old
                end

            end

        end

 
    end

    if model_spec.symmetric_W
        symmetrize!(W)
    end

    inds = get_augmented_submatrix_indices(model.augmented_weights, 1, 0)
    model.augmented_weights[inds,inds] = W

end

function sample_W_full(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       Z,#sparse or full 2d array
                       observed_effects::Array{Float64,2})

    YY = data.Ytrain
    (K,K) = size(model.weights)
    (N,N) = size(YY[1])
    W = copy(model.weights)

    num_W_sweeps = 3

    if model_spec.symmetric_W
        A = ones(K,K)
        triu_inds = find(triu(A))
        vectorize = x -> x[triu_inds]
        devectorize = x -> (A = zeros(K,K); A[triu_inds] = x; symmetrize(A)) 
    elseif model_spec.diagonal_W
        vectorize = x -> diag(x)
        devectorize = x -> diagm(x)
    else
        vectorize = x -> x[:]
        devectorize = x -> reshape(x, (K,K))
    end

    temp_model = copy(model) 

    function W_logpdf_vectorized(x::Vector{Float64})
        W_matrix = devectorize(x)
        temp_model.weights = W_matrix
        W_full_logpdf(temp_model, model_spec, data, Z, observed_effects)
    end
    function W_logpdf_gradient_vectorized(x::Vector{Float64})
        W_matrix = devectorize(x)
        temp_model.weights = W_matrix
        grad = W_full_logpdf_gradient(temp_model, model_spec, data, Z, observed_effects)
        vectorize(grad)
    end

    hmc_opts = @options L=4 stepsize=0.02

    W = vectorize(W) 
    for iter = 1:num_W_sweeps
        println("W sampling sweep ", iter,"/",num_W_sweeps)
        
        W = hmc_sampler(W, W_logpdf_vectorized, W_logpdf_gradient_vectorized, hmc_opts)
    end

    W = devectorize(W)
    model.weights = W

    inds = get_augmented_submatrix_indices(model.augmented_weights, 1, 0)
    model.augmented_weights[inds,inds] = W
end


function sample_Z(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState,
                  latent_effects::Array{Float64,2},
                  observed_effects::Array{Float64,2},
                  model_logprob::Float64)

    println("Sample Z")
    YY = data.Ytrain
    (N,N) = size(YY[1])
    tree = model.tree
    num_W_sweeps = 3
    num_W_slice_steps = 1
    # Sample new features from root to leaf
    root = FindRoot(tree, 1) 
    leaf_to_root = GetLeafToRootOrdering(tree, root.index)
    root_to_leaf = reverse(leaf_to_root)


    mu = ones(2N - 1)
    for i = root_to_leaf
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)
            cur_mu_prop = self_direction == 1 ? (1-parent.rho) : parent.rho
            cur_mu_prop *= parent.rhot
            mu[i] = mu[parent.index]*cur_mu_prop
        end
    end


    Z = ConstructZ(tree)
    (K, K) = size(model.weights)
    #Array(Array{Int64, 2}, (K,K))

    old_relevant_pairs = compute_all_relevant_pairs(model_spec, data, K, Z)

    current_model_logprob = model_logprob

    for rtl_index = 1:length(root_to_leaf)
        tree = model.tree #make sure after copies to model we have the right tree ref
        node_index = root_to_leaf[rtl_index]
        
        # root node cannot have features above it
        if node_index == root.index
            continue
        end
        parent_index = tree.nodes[node_index].parent.index

        if mod(rtl_index,ceil(length(root_to_leaf)/10)) == 0
            percent = ceil(rtl_index/ceil(length(root_to_leaf)/10))*10
            println(" ",percent , "% ")
        end
#        if model_spec.symmetric_W
#            println("rtl_index: ", rtl_index)
#            (K0,K0) = size(model.weights)
#            if K0 > 1
#                println(old_relevant_pairs[1,K0])
#                assert( length(old_relevant_pairs[1,K0]) == 0)
#            end
#            println("rtl_index: ", rtl_index)
#
#        end
#
        if rand() > model_spec.Z_sample_branch_prob
            continue
        end

        old_model_logprob = current_model_logprob

        if model_spec.verbose
            println("vardim sampling node_index: ", node_index)
        end
        u = tree.nodes[node_index].state

        tU = [tree.nodes[i].state for i = 1:2N-1]

        if !all(model.augmented_weights.num_active_features .== tU)
            fail_inds = find( model.augmented_weights.num_active_features .!= tU)
            println(fail_inds)
            println(model.augmented_weights.num_active_features[fail_inds])
            println(tU[fail_inds])

            println(model.augmented_weights.num_active_features[end-20:end])
            println(tU[end-20:end])
            println(sum(model.augmented_weights.num_active_features))
            println(sum(tU))
        end
        assert(all(model.augmented_weights.num_active_features .== tU))

        effective_lambda = model.lambda * (mu[parent_index]^model.gamma - mu[node_index]^model.gamma)

        (K, K) = size(model.weights)
        
        inds = get_augmented_submatrix_indices(model.augmented_weights, 1, 0)
        augW = model.augmented_weights[inds,inds] 

#        println("numfeatures: ", sum(model.augmented_weights.num_active_features))
#        println("size(augW): ", size(augW))
#        println("K: ", K)
        if !all(model.weights[1:K,1:K] .== augW)
            println("model_inds: ", inds)
            println(model.weights[1:K,1:K])
            println(augW)
            println(model.augmented_weights.matrix)
        end

        assert(all(model.weights[1:K,1:K] .== augW))


        new_model = copy(model)
        L = u + randmult(model_spec.rrj_jump_probabilities) - 2

        # if L < 0, then we are left with a single, discrete component in s_L, so 
        # no need to sample, we remain in the same state u = 0
        if L < 0
            continue
        end

        if model_spec.debug
            println("(L,u,ind): ", (L,u,rtl_index))
            Z = ConstructZ(model.tree)
            W = copy(model.weights)
            symmetrize!(W)
            verify_latent_effects = Z*W*Z' 
            badinds = findn(abs(verify_latent_effects - latent_effects) .> 10.0^-5)

            println(badinds)
            println(latent_effects[badinds...])
            println(verify_latent_effects[badinds...])
            println(model.weights)
            println(W)
            assert(max(abs(verify_latent_effects - latent_effects)) < 10.0^-5)

        end

        oldK = K
        K = K + L - u

        new_model.tree.nodes[node_index].state = L

        UU = [new_model.tree.nodes[x].state for x = 1:length(new_model.tree.nodes)]
        assert(sum(UU) == K)

        old_W_index_pointers = weight_index_pointers(model.tree)
        new_W_index_pointers = weight_index_pointers(new_model.tree)
        start_index = new_W_index_pointers[node_index]
        end_index = start_index + L - 1
        new_k  = end_index + 1


        assert(old_W_index_pointers[node_index] == new_W_index_pointers[node_index])
        old_latent_effects = copy(latent_effects)

        aug_k = 0
        # pick one feature to remove to be the auxiliary for
        # s_{u-1}'s k term

        st_ind = old_W_index_pointers[node_index]
        end_ind = st_ind + u - 1

        adj_relevant_pairs = compute_new_relevant_pairs(model,model_spec,data,
                             old_relevant_pairs, node_index, st_ind, end_ind)

        if L == u-1
            new_relevant_pairs = copy(old_relevant_pairs)

            aug_u = rand(1:u)
            deactivate_feature(new_model.augmented_weights,
                               node_index,
                               aug_u)

       
            aug_k = old_W_index_pointers[node_index] + aug_u - 1
            aug_new_k = old_W_index_pointers[node_index] + u - 1

            permutation = [1:oldK]
            permutation[aug_k] = aug_new_k
            permutation[aug_new_k] = aug_k
            
            permute_rows_and_cols!(new_relevant_pairs, permutation)
            new_relevant_pairs = delete_row_and_col(new_relevant_pairs, aug_new_k)

            latent_effects = adjust_latent_effects(model, model_spec,
                                 latent_effects, old_relevant_pairs, u, L, aug_k, 
                                 0)

            aug_k = new_W_index_pointers[node_index] + aug_u - 1
        elseif L == u
            new_relevant_pairs = copy(old_relevant_pairs)
        elseif L == u+1
            new_relevant_pairs = copy(adj_relevant_pairs)  
        
            activate_feature(new_model.augmented_weights,
                             node_index) 

            new_model_inds = 
                get_augmented_submatrix_indices(new_model.augmented_weights, node_index, 0)
            new_model.weights = new_model.augmented_weights[new_model_inds, new_model_inds]
            aug_k = end_ind + 1 
            latent_effects = adjust_latent_effects(new_model, model_spec, latent_effects,
                                 new_relevant_pairs, u, L, aug_k, aug_k)
           
        end


        if model_spec.debug

            println("starting model_logprob: ", current_model_logprob)
            tmodel = copy(model)
            tmodel.tree.nodes[node_index].state = L
            Z = ConstructZ(tmodel.tree)

            verify_relevant_pairs = compute_all_relevant_pairs(model_spec, data, K, Z)

            assert( verify_relevant_pairs == new_relevant_pairs)


            WWinds = 
                get_augmented_submatrix_indices(model.augmented_weights, node_index, 1)
            WW = model.augmented_weights[WWinds,WWinds]
            tmodel = copy(model)
            tmodel.weights = WW
            verify_latent_effects = compute_latent_effects(tmodel, model_spec, u, L,
                                        aug_k, node_index)

            permutation = sortperm( abs(latent_effects - verify_latent_effects)[:])
            rev = reverse(permutation)
            maxinds = rev[1:10]
          
#            println("Latent Effects:")
#            println(ind2sub(size(latent_effects), maxinds) )
#            println(latent_effects[maxinds])
#            println(verify_latent_effects[maxinds])
#            println(latent_effects[maxinds] - verify_latent_effects[maxinds])
#            println(latent_effects[1:10] - verify_latent_effects[1:10])

            println("L,u,aug_k: ", (L,u,aug_k)) 
            assert(max(abs(latent_effects - verify_latent_effects)) < 10.0^-5)
        end

        new_model_inds = 
            get_augmented_submatrix_indices(new_model.augmented_weights, node_index, 1)
        new_model.weights = new_model.augmented_weights[new_model_inds, new_model_inds] 

        new_relevant_pairs = compute_new_relevant_pairs(new_model,model_spec,data,
                                 new_relevant_pairs, node_index, start_index, end_index)


        if model_spec.debug
            tmodel = copy(model)
            tmodel.tree.nodes[node_index].state = L+1
            Z = ConstructZ(tmodel.tree)

            verify_relevant_pairs = compute_all_relevant_pairs(model_spec, data, K+1, Z)
            assert( verify_relevant_pairs == new_relevant_pairs)
        end

        # adjust model probability
        adjust_terms = adjust_model_logprob(model, model_spec, data, adj_relevant_pairs,
                           old_latent_effects, latent_effects, observed_effects, aug_k,
                           u, L, old_W_index_pointers, node_index, effective_lambda)

        assert( L == u ? sum(adjust_terms[1:2]) == 0.0 : true)

        current_model_logprob += sum(adjust_terms[1:2])

        W = new_model.weights


#        println("L,u: ", (L,u))
#        println("end_index: ", end_index)
#        println("size(new_relevant_pairs): ", size(new_relevant_pairs))
#        UU = [new_model.tree.nodes[x].state for x = 1:length(new_model.tree.nodes)]
#        println("sumU: ",sum(UU))




        component_latent_effects = adjust_component_latent_effects(new_model, model_spec, 
                                       latent_effects, new_relevant_pairs, L, 
                                       start_index, end_index)

        k_is_auxiliary = Array(Array{Bool,1},0)
 
        # Construct w_is_auxiliary for each component
        w_is_auxiliary = Array(Array{Bool,2},0)

        num_local_mutations = [ x <= L ? L-1 : x == L+1 ? L : L+1 for x = 1:L+2] # [u-1, u-1, u-1, ..., u-1, u, u+1]

        if model_spec.debug
            verify_latent_effects = compute_component_latent_effects(new_model,
                                        model_spec, L, start_index, end_index, node_index)

            le_diffs = verify_latent_effects .- component_latent_effects

            maxabss = [ max(abs(le_diffs[x])) for x = 1:length(le_diffs)]
            println(maxabss )
            println(maxabss .< 10.0^-5)


            tmodel = copy(new_model)
            tmodel.tree.nodes[node_index].state = L+1
            ZZ = ConstructZ(tmodel.tree)
            bad_indices = find(maxabss .> 10.0^-5)
            if length(bad_indices) > 0

                bad_ind = bad_indices[1]
                println("bad_ind: ", bad_ind)
                I,J = findn(abs(verify_latent_effects[bad_ind] - component_latent_effects[bad_ind]) .> 10.0^-5)


                println(size(ZZ))
                println(I[1])
                println(J[1])
          
                println(ZZ[I[1],:])
                println(ZZ[J[1],:]) 

                println(new_model.weights)
                maxi,maxj = findn(maxabss[bad_ind] .== abs(le_diffs[bad_ind]))
                println("original: ", latent_effects[maxi[1],maxj[1]])
                println("verified: ", verify_latent_effects[bad_ind][maxi[1],maxj[1]])
                println("component: ", component_latent_effects[bad_ind][maxi[1],maxj[1]])


                k_common = find(ZZ[I[1],:] .* ZZ[J[1],:] .> 0)
                println("k_common: ", k_common)
                println(W[k_common, k_common])

                k1 = k_common[1]
                k2 = k_common[2]
                II = find( new_relevant_pairs[k1,k2][1,:] .== I[1])
                println(II)
                println(new_relevant_pairs[k1,k2][:,II])
            end

            assert(all(maxabss .< 10.0^-5))

        end
 
        for u_ind = 1:L
            k_ind = start_index + u_ind - 1
            push!(w_is_auxiliary,[x == k_ind  || y == k_ind ||
                                 x == new_k || y == new_k
                                 for x = 1:K+1, y = 1:K+1] )

            push!(k_is_auxiliary, [x == k_ind || x == new_k for x = 1:K+1])
            
        end

        push!(w_is_auxiliary,[x == new_k || y == new_k for x = 1:K+1, y = 1:K+1] )
        push!(k_is_auxiliary,[x == new_k for x = 1:K+1] )

        push!(w_is_auxiliary, zeros(Bool, (K+1,K+1)))
        push!(k_is_auxiliary, zeros(Bool, K+1))

        # "constant" terms needed so that the sum of the different mixture components
        # is done with each scaled correctly when performing local pdf computations
        # ie conditionals are not enough when the total mixture is in the joint space
        const_terms = compute_constant_terms(new_model, model_spec, data, 
                          new_relevant_pairs, component_latent_effects,
                          observed_effects, num_local_mutations, new_W_index_pointers,
                          node_index, effective_lambda, current_model_logprob)


        if model_spec.debug
            const_terms2 = vardim_splits(new_model, model_spec, data, component_latent_effects,
                              observed_effects, num_local_mutations, new_W_index_pointers, node_index,
                              effective_lambda, w_is_auxiliary)

            const_terms3 = const_terms + 
                vardim_multiplier_terms(model_spec, num_local_mutations)

            println("adjusted const terms:")
            println(const_terms3)
            println("vardim_splits const terms:")
            println(const_terms2)
            if abs(max(abs(const_terms2 - const_terms3))) > 0.1
                println("SLDKJF")
                println(current_model_logprob)
                println(exp_normalize(const_terms3))
                println(exp_normalize(const_terms2))
                println(const_terms3 - max(const_terms3))
                println(const_terms2 - max(const_terms2))
            end
            curind = length(const_terms)-1
            unscaled = const_terms
                #const_terms - vardim_multiplier_terms(model_spec, num_local_mutations)

            verify_state_logprob = compute_unaugmented_prob(new_model, model_spec, 
                                          new_relevant_pairs, component_latent_effects,
                                          observed_effects, num_local_mutations,
                                          new_W_index_pointers, node_index,
                                          effective_lambda, curind,
                                          unscaled[curind])
            if abs(verify_state_logprob - current_model_logprob) > 10.0^-5
                println("verify prob is off")
                println(verify_state_logprob)
                println(current_model_logprob)
            end
        end


        if model_spec.debug && model_spec.verbose
            println("W[start,start]: ", new_model.weights[start_index,start_index])
            println("slice sample")
        end

        old_const_terms = copy(const_terms)

        for iter = 1:num_W_sweeps
            k1_range = start_index:end_index+1 #model_spec.diagonal_W ? new_k : start_index:end_index+1 #(1:K+1)
            for k1 = k1_range

                k2_range = model_spec.diagonal_W ? k1 :
                           model_spec.symmetric_W ? (start_index:k1) :
                           (start_index:end_index+1)

                for k2 = k2_range

                    #println("slice sample prep, start, end:", start_index, " ", end_index)
                    w_is_auxiliary_k1k2 = zeros(Bool, length(num_local_mutations))

                    for uuu = 1:length(w_is_auxiliary_k1k2)
                        w_is_auxiliary_k1k2[uuu] = w_is_auxiliary[uuu][k1,k2]
                    end

                    w_cur = W[k1,k2]
                    w_old = w_cur

#                    println("size(new_relevant_pairs): ", size(new_relevant_pairs))
#                    println("k1: ", k1)
#                    println("k2: ", k2)
                    delta_const_terms = vardim_local_splits(new_model, model_spec, data, 
                        new_relevant_pairs[k1,k2], component_latent_effects,
                        observed_effects, num_local_mutations, new_W_index_pointers,
                        node_index, effective_lambda, zeros(length(num_local_mutations)),
                        w_old, w_old, w_is_auxiliary_k1k2) 

                    const_terms -= delta_const_terms

                    g = x -> vardim_local_sum(new_model, model_spec, data,
                                 new_relevant_pairs[k1,k2],
                                 component_latent_effects, observed_effects,
                                 num_local_mutations, new_W_index_pointers, node_index,
                                 effective_lambda, const_terms, w_old, x,
                                 w_is_auxiliary_k1k2) 
                    gx0 = g(w_cur)
                   
                    for slice_iter = 1:num_W_slice_steps
                        (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                    end
                    W[k1,k2] = w_cur

#                    if model_spec.debug
#                        delta_const_terms = vardim_local_splits(new_model, data, 
#                            new_relevant_pairs[k1,k2], component_latent_effects,
#                            observed_effects, num_local_mutations, W_index_pointers,
#                            node_index, effective_lambda, zeros(length(num_local_mutations)),
#                            w_old, w_cur, w_is_auxiliary_k1k2)
#                        println("w_old, w_new, prob: ", w_old,", ", w_cur, ", ", g(w_cur))
#                        println(delta_const_terms)
#                        println(const_terms)
#                    end

                    delta_const_terms = vardim_local_splits(new_model, model_spec, data, 
                        new_relevant_pairs[k1,k2], component_latent_effects,
                        observed_effects, num_local_mutations, new_W_index_pointers,
                        node_index, effective_lambda, zeros(length(num_local_mutations)),
                        w_old, w_cur, w_is_auxiliary_k1k2) 

                    const_terms += delta_const_terms
        
#                    println("delta_const_terms: ", delta_const_terms)
#                    println("num_local_mutations: ", num_local_mutations)
#                    println("k1,k2: ", (k1,k2))
#                    println(const_terms)
#                    println(old_const_terms)
#                    println(const_terms - old_const_terms)
#                    assert(const_terms[1] == old_const_terms[1] || num_local_mutations[1] > 0)


                    #println("adjust precomputations")
                    rpairs = new_relevant_pairs[k1,k2]
                    u_ind_range = find(!w_is_auxiliary_k1k2)
                    for p = 1:size(rpairs)[2]
                        i = rpairs[1,p]
                        j = rpairs[2,p]

                        latent_effects[i,j] += w_cur - w_old


                        for u_ind = u_ind_range 
                            component_latent_effects[u_ind][i,j] += w_cur - w_old
                        end

                    end


#                    unscaled_const_terms = 
#                        const_terms - vardim_multiplier_terms(model_spec, num_local_mutations)
#
#                    println(unscaled_const_terms) 
                end
            end
        end

        if model_spec.symmetric_W
            symmetrize!(W)
        end

        logprobs = const_terms

        scaled_logprobs = 
            logprobs + vardim_multiplier_terms(model_spec, num_local_mutations)

        if model_spec.debug
            logprobs2 = vardim_splits(new_model, model_spec, data, component_latent_effects,
                                     observed_effects, num_local_mutations, new_W_index_pointers, node_index,
                                     effective_lambda, w_is_auxiliary)

            if norm(exp_normalize(scaled_logprobs) - exp_normalize(logprobs2)) > 10.0^-5
                println("BLAH")
                println(scaled_logprobs)
                println(logprobs2)
            end

            println("const_terms: ")
            println(const_terms)
            println("old_const_terms: ")
            println(old_const_terms)
        end

        sampled_component = randmult(exp_normalize(scaled_logprobs)) 

        new_u = num_local_mutations[sampled_component]
        new_K = K + new_u - L
        new_valid_indices = find(!w_is_auxiliary[sampled_component])


#        println("new_u: ", new_u)
#        println(w_is_auxiliary[sampled_component])
#        println(new_valid_indices)
#        println("new_K: ",new_K)
#        println("K: ",K)
#        println("L: ",L)
        new_W = reshape(copy(W[new_valid_indices]), (new_K, new_K)) 

        oldW = copy(model.weights)
        augW = copy(W)


        proposed_model = copy(new_model)
        proposed_model.tree.nodes[node_index].state = new_u
        proposed_model.weights = copy(new_W)


        proposed_model_inds = 
            get_augmented_submatrix_indices(proposed_model.augmented_weights,
                                            node_index,
                                            1)
#        println(new_valid_indices)
#        println("sampled_component: ", sampled_component)
#        println(w_is_auxiliary[sampled_component])
#        println("proposed_inds: ", proposed_model_inds)
#        println(W)
#        println(new_W)

        proposed_model.augmented_weights[proposed_model_inds, proposed_model_inds] =
            augW 

        

        if new_K < K
            deactivate_feature(proposed_model.augmented_weights,
                               node_index,
                               sampled_component)
        elseif new_K > K
            activate_feature(proposed_model.augmented_weights,
                             node_index)
        end

        next_inds = get_augmented_submatrix_indices(proposed_model.augmented_weights,
                                                    node_index,
                                                    0)
        if model_spec.debug
            assert(proposed_model.augmented_weights[next_inds,next_inds] == new_W)
        end
#        permutation = [new_valid_columns, invalid_columns]
#        permute_rows_and_cols!(augW, permutation)
#        permute_rows_and_cols!(proposed_model.augmented_weights, permutation)
#        if augmented_K < K+1
#            proposed_model.augmented_weights = [zeros(2N-1) for k1 = 1:K+1, k2 = 1:K+1]
#        else
#            proposed_model.augmented_weights[1:K+1,1:K+1] = augW
#        end
#
#
#        for k1 = 1:new_K
#            for k2 = 1:new_K
#                proposed_model.augmented_weights[k1,k2][:] = proposed_model.weights[k1,k2]
#            end
#        end



        proposed_state_logprob = compute_unaugmented_prob(new_model, model_spec, 
                                      new_relevant_pairs, component_latent_effects,
                                      observed_effects, num_local_mutations,
                                      new_W_index_pointers, node_index,
                                      effective_lambda, sampled_component,
                                      logprobs[sampled_component])

        inds = get_augmented_submatrix_indices(proposed_model.augmented_weights, node_index, 0)
        augW = proposed_model.augmented_weights[inds,inds]

 
#        println("inds: ", inds)
#        println(inds)
#        println("prop_inds: ", proposed_model_inds)
#        println(proposed_model_inds)
#        println(W)
#        println(w_is_auxiliary[sampled_component])
#        println(augW)
#        println(proposed_model.weights)
        assert(all(proposed_model.weights .== augW))

        # acceptance prob 1
        model.weights = copy(proposed_model.weights)
        model.tree = copy(proposed_model.tree)
        model.augmented_weights = copy(proposed_model.augmented_weights)


#        println("size(new_relevant_pairs): ", size(new_relevant_pairs))
#        println("new_valid_indices: ", new_valid_indices)
        old_relevant_pairs = reshape(new_relevant_pairs[new_valid_indices], (new_K, new_K))
        latent_effects = component_latent_effects[sampled_component]
        current_model_logprob = proposed_state_logprob 

        UU = [model.tree.nodes[x].state for x = 1:length(model.tree.nodes)]
        assert(new_K == sum(UU))




#        println("new_u,weight_prob,new_weight_prob ", new_u, ",", sum(normal_logpdf(oldW, model.w_sigma )), ",", sum(normal_logpdf(new_W, model.w_sigma )))

        if model_spec.verbose
            println("oldK,newK,L,u,eff_lambda,Delta_K: ", oldK,",",new_K, ",",L,",",u,",", effective_lambda, ",",new_K - oldK)
            println(logprobs)
            println(scaled_logprobs)
            println(num_local_mutations)
            println("model_logprob: ", proposed_state_logprob)
        end


        if model_spec.debug 

            Z = ConstructZ(model.tree)

            full_latent_effect = Z * model.weights * Z'

            println("maxdiff latent effect: ", max(abs(full_latent_effect-latent_effects)))

            full_probs = zeros(length(logprobs))
            for cmp = 1:length(logprobs)

                new_uu = num_local_mutations[cmp]
                new_K = K + new_uu - L
                new_valid_indices = find(!w_is_auxiliary[cmp])
                new_W = reshape(copy(W[new_valid_indices]), (new_K, new_K)) 

                new_model.weights = copy(new_W)
                new_model.tree.nodes[node_index].state = new_uu

                tree_prior = prior(new_model,model_spec) 
                tree_LL, test_LL = likelihood(new_model, model_spec, data, [1:N])
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

        if model_spec.debug && new_u - L > 0
            tree_prior = prior(model,model_spec) 
            tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
        end



        if model_spec.debug && new_u - L > 0
            if abs(current_model_logprob - (tree_prior + tree_LL)) > 10.0^-5
                println((start_index,end_index,new_k))
                println(current_model_logprob)
                println(tree_prior + tree_LL)
                println(new_relevant_pairs[new_k,new_k])
                println(new_relevant_pairs[start_index,start_index])
                 
                assert(length(new_relevant_pairs[new_k,new_k]) > 0)
            end
        end

    end
end

function sample_rho_rhot(model::ModelState,
                         model_spec::ModelSpecification,
                         slice_iterations::Int)

    println("Sample rho, rhot")
    tree = model.tree
    gam = model.gamma
    lambda = model.lambda

    N::Int = (length(tree.nodes)+1)/2

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)
    T = ones(2N-1)
    S_f = zeros(2N-1)
    K = zeros(2N-1)

    for i = reverse(indices)
        if i > N
            cur = tree.nodes[i]
            l = tree.nodes[i].children[1].index
            r = tree.nodes[i].children[2].index
            if i != root.index
                left_prob = tree.nodes[i].rhot*tree.nodes[i].rho
                right_prob = tree.nodes[i].rhot*(1.0-tree.nodes[i].rho)

                T[l] = T[i] * left_prob^gam
                T[r] = T[i] * right_prob^gam
            end
        end
    end

    for i = indices
        if i > N
            left_child = tree.nodes[i].children[1]
            right_child = tree.nodes[i].children[2]
            K[i] = K[left_child.index] + K[right_child.index] +
                   left_child.state + right_child.state 
       
            left_prob = (1-tree.nodes[i].rhot*tree.nodes[i].rho)^gam
            right_prob = (1-tree.nodes[i].rhot*(1.0-tree.nodes[i].rho))^gam
 
            S_f[i] = S_f[left_child.index] + S_f[right_child.index] + 
                     left_prob*T[left_child.index] + right_prob*T[right_child.index] 
        end
    end

    T_mult = ones(2N-1)

    for i = reverse(indices)
        if i > N
            cur = tree.nodes[i]
            left_child = tree.nodes[i].children[1]
            right_child = tree.nodes[i].children[2]
            l = left_child.index
            r = right_child.index

            k_l = left_child.state
            k_r = right_child.state
            N_l = left_child.num_leaves
            N_r = right_child.num_leaves
             
            rho = cur.rho
            rhot = cur.rhot
            nut_l = cur.rhot
            nut_r = cur.rhot
            nu_l = cur.rho
            nu_r = 1-cur.rho

            T_l = T[l] * T_mult[i]
            T_r = T[r] * T_mult[i]

            S_l = S_f[l] * T_mult[i] / (nut_l * nu_l)^gam
            S_r = S_f[r] * T_mult[i] / (nut_r * nu_r)^gam 
     
            f = x -> rho_logpdf(x, gam, lambda, nut_l, nut_r, k_l, k_r, 
                                K[l], K[r], T[l], T[r], S_l, S_r, N_l, N_r) 
            (rho, f_rho) = slice_sampler(rho, f, 0.1, 10, 0.0, 1.0)


            p_s = 2/(N_l+N_r+1)
            f = x -> logsumexp( rhot_splits(x, p_s, gam, lambda, nu_l, nu_r, k_l, k_r,
                                      K[l], K[r], T[l], T[r], S_l, S_r))
            (rhot_u, f_rhot) = slice_sampler(rhot, f, 0.1, 10, 0.0, 1.0)

            f_vals = rhot_splits(rhot_u, p_s, gam, lambda, nu_l, nu_r, k_l, k_r,
                                 K[l], K[r], T[l], T[r], S_l, S_r)

            f_ind = rand(Categorical(exp_normalize(f_vals)))
            rhot = f_ind == 1 ? 1.0 : rhot_u

            # Update for T and S
            l_new = rhot*rho
            r_new = rhot*(1-rho)
            l_old = cur.rhot*cur.rho
            r_old = cur.rhot*(1-cur.rho) 

            T_mult[l] = T_mult[i] * (l_new/l_old)^gam
            T_mult[r] = T_mult[i] * (r_new/r_old)^gam


            cur.rho = rho
            cur.rhot = rhot
        end
    end

end



function sample_psi(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState,
                    latent_effects::Array{Float64,2},
                    observed_effects::Array{Float64,2})

    tree = model.tree
    N::Int = (length(tree.nodes)+1)/2

    println("psi: ")
    for prune_index = 1:N # only prune leaves, it's much faster

        if mod(prune_index,ceil(N/10)) == 0
            percent = ceil(prune_index/ceil(N/10))*10
            println(" ",percent , "% ")
        end

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

        if parent.children[1].index == prune_index
            original_sibling = parent.children[2]
        else
            original_sibling = parent.children[1]
        end

        prune_tree!(model, prune_index)
        gp = grandparent.index




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

        if true
            psi_p = prior_tree(tree, prune_index)
        end 
        (priors, pstates) = psi_infsites_logpdf(model, data, prune_index, path)
        (likelihoods, lstates) = psi_likelihood_logpdf(model, model_spec, data, prune_index, path)

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)

        if any(isnan(probs))
            nan_ind = find(isnan(probs))[1]
            maxprior = max(priors)
            maxll = max(likelihoods)
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
                push!(A, i)
            end 

            push!(graft_indices, pstates[i][1])
            push!(new_graftpoint_features, pstates[i][2])
            push!(new_parent_features, pstates[i][3])
        end



        assert(pstates[state_index] == lstates[state_index])

        (graft_index, graftpoint_features, parent_features) = pstates[state_index]

        graft_tree!(model, prune_index, graft_index, parent_features, graftpoint_features)


        if true #model_spec.debug 
            println("Sampling Prune Index: ", prune_index, " Num Leaves: ", length(GetLeaves(tree, grandparent.index)))
            println("Num Leaves pruned: ", length(GetLeaves(tree, prune_index)), " Num leaves remaining: ", length(GetLeaves(tree, gp)) )
            println("original_index,insert_index,parent,root: ", original_sibling.index, ",", pstates[state_index][1], ",", parent.index, ",", root.index)
#            println("logprobs: ", logprobs)
#            println("graft indices: ", graft_indices)
#            println("new parent features: ", new_parent_features)
#            println("new graftpoint features: ", new_graftpoint_features)
            println("start_probs, end_prob: ",logprobs[A]  ,", ", logprobs[state_index])
            println("start_prior, start_LL: ",priors[A], ", ", likelihoods[A])
            println("end_prior, end_LL: ", priors[state_index], ", ", likelihoods[state_index])

            psi_before = psi_p[original_sibling.index]
            psi_after = psi_p[graft_index]

            println("psi_diff: $(psi_after-psi_before)")

            tree_prior = prior(model, model_spec) 
            tree_LL, test_LL = likelihood(model, model_spec, data, [1:N])
            println("tree probability, prior, LL, testLL: ", tree_prior + tree_LL, ",", tree_prior, ",",tree_LL, ",", test_LL)
        end
    end
end

function construct_effects(model::ModelState,
                           model_spec::ModelSpecification,
                           data::DataState,
                           Z) #sparse or full binary array

    YY = data.Ytrain
    W = model.weights
    (N,N) = size(YY[1])
    #latent_effects = zeros(Float64, (N,N))
    observed_effects = zeros(Float64, (N,N))
    for i = 1:N
        for j = 1:N
            #latent_effects[i,j] = (Z[i,:] * W * Z[j,:]')[1]
            observed_effects[i,j] = compute_observed_effects(model, model_spec, data, i, j)
        end
    end
    latent_effects = Z*W*Z'
    (latent_effects, observed_effects)

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
