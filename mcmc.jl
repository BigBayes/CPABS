
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
function mcmc_sweep(tree::Tree{Float64},
                    full_move_prob::Float64,
                    model::ModelParameters)
    N = (length(tree.nodes)+1)/2
    int_exp_elapsed = 0.0
    sample_elapsed = 0.0
    for prune_index = 1:2N-2
        parent = tree.nodes[prune_index].parent
        grandparent = tree.nodes[prune_index].parent.parent
        if grandparent == Nil{Float64}()
            continue
        end
        PruneIndexFromTree!(tree,prune_index)
        gp = grandparent.index

        u = rand()
        if u > full_move_prob
            descendant_ind = GetRandomDescendant(tree,gp)
            prior_probs = prior_path(tree,prune_index,descendant_ind)

            path = GetPath(tree,descendant_ind) 
        else
            path = tree.nodes
        end        

 
        funcs = cell(length(path))
        Z = zeros(length(path))
        L = -Inf*ones(length(path))
        A = -ones(length(path))
        B = -ones(length(path))
        for i = 1:length(path)-1 #Don't add above the root
            f = continuous_pdf(tree,prune_index,path[i].index,prior_probs[i],model)
            a = max(path[i].coalescent_time,
                    tree.nodes[prune_index].coalescent_time)
            b = path[i].parent.coalescent_time
            if a < b
                funcs[i] = f
                A[i] = a
                B[i] = b
                s = time()
                Z[i],L[i] = int_exp(f,a,b)
                int_exp_elapsed += time()-s
            end
        end
        Lold = L
        Zold = Z
        L = L - max(L)
        Z = Z .* exp(L)
        path_index = rand_multinomial(Z)
#            if path_index == 0
#                global F = funcs
#                global T = tree
#                global P = path
#                global A = A
#                global B = B
#                global Lold = Lold
#                global Zold = Zold
#            end 
        a = max(path[path_index].coalescent_time,
                tree.nodes[prune_index].coalescent_time)
        b = path[path_index].parent.coalescent_time
        s = time()
        t = sample_time(funcs[path_index], a, b, 100, 0.2) 
        sample_elapsed += time()-s

        InsertIndexIntoTree!(tree, prune_index, path[path_index].index)
        parent.coalescent_time = t

        neighbor_locations = [parent.parent.location',
                              parent.children[1].location',
                              parent.children[2].location']
        neighbor_variances = model.k*[parent.parent.coalescent_time - t,
                              t - parent.children[1].coalescent_time,
                              t - parent.children[2].coalescent_time]

        parent.location = rand_normprod(neighbor_locations,
                                        neighbor_variances)[:]
        
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
