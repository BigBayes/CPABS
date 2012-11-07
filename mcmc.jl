load("model.jl")
load("tree.jl")
load("probabiltiy_util.jl")

function mcmc(Y::Array{Int64,2},
              X_r::Array{Float64,3},
              X_p::Array{Float64,2},
              X_c::Array{Float64,2},
              K_init::Int64,
              lambda::Float64,
              gamma::Float64,
              iterations::Int64)
    (N1, N2, pairwise_p) = size(X_r)
    assert(N1 == N2)
    N = N1
    (N1, parent_p) = size(X_p)
    (N2, children_p) = size(X_c)
    assert(N1 == N)
    assert(N2 == N)

    U = zeros(Int64,2N-1)
    for i = 1:2N-1
        U[i] = randi(3)-1
    end

    W = randn(sum(U),sum(U))
    tree = Tree(U)

    model = ModelState(lambda,gamma,tree,W)

    for iter = 1:iterations
        prior = tree_prior(tree)
        branch_prob = tree_branch_prob(tree, model)
        data_prob = tree_likelihood(tree, model)
        println("Tree probability: ",(prior + branch_prob + data_prob,
                                      prior, branch_prob, data_prob))
        mcmc_sweep(tree, 0.0, model)
    end

end


#@profile begin
function mcmc_sweep(model::ModelParameters,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2})

    tree = model.tree
    N = (length(tree.nodes)+1)/2
    int_exp_elapsed = 0.0
    sample_elapsed = 0.0
    for prune_index = 1:2N-1
        parent = tree.nodes[prune_index].parent
        if parent == Nil()
            continue
        end

        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil()
            continue
        end

        prune_tree!(tree, prune_index)
        gp = grandparent.index


        (priors, pstates) = infsites_logpdf(model, prune_index)
        (likelihoods, lstates) = psi_likelihood_logpdf(model, prune_index, Y, X_r, X_p, X_c)

        logprobs = priors + likelihoods
        probs = exp_normalize(logprobs)

        state_index = randmult(probs)

        (graft_index, graftpoint_features, parent_features) = pstates[state_index]

        graft_tree!(model, prune_index, graft_index, parent_features, graftpoint_features)
    end
    println(int_exp_elapsed, sample_elapsed)
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
