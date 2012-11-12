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
        U[i] = randi(3)-1
    end

    W = randn(sum(U),sum(U))
    tree = Tree(U)

    model = ModelState(lambda,gamma,1.0,tree,W,[0.0],[0.0],[0.0],[0.0],[0.0],0.0)

    for iter = 1:iterations
        println("Iteration: ", iter)
        tree_prob = full_pdf(model, model_spec, Y, X_r, X_p, X_c)
        println("Tree probability: ", tree_prob)
        mcmc_sweep(model, model_spec, Y, X_r, X_p, X_c)
    end
end


@profile begin
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
    for prune_index = 1:10

        println("Sampling Prune Index ", prune_index)
        parent = tree.nodes[prune_index].parent
        if parent == Nil()
            continue
        end

        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil()
            continue
        end

        prune_tree!(model, prune_index)
        gp = grandparent.index


        (priors, pstates) = psi_infsites_logpdf(model, prune_index)
        (likelihoods, lstates) = psi_likelihood_logpdf(model, prune_index, Y, X_r, X_p, X_c)

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
   
        for k1 = 1:3
        #for k1 = 1:K
            println("W sampling k1: ", k1)
            for k2 = 1:3
            #for k2 = 1:K
                relevant_pairs = compute_relevant_pairs(Z,k1,k2)
                g = x -> W_local_logpdf(model, Y, relevant_pairs, latent_effects, observed_effects, k1, k2, x)

                w_cur = W[k1,k2]
                gx0 = g(w_cur)
                for slice_iter = 1:num_W_slice_steps
                    (w_cur, gx0) = slice_sampler(w_cur, g, 1.0, 10, -Inf, Inf, gx0)
                end

                W[k1,k2] = w_cur
            end

        end 
    end

    W_time = time() - W_time

    println("MCMC Timings (psi, W) = ", (psi_time, W_time))
end

end #profile

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
