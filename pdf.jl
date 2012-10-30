#Global pdf
function full_pdf(model::ModelState,
                  Y::Array{Int64,2},
                  X_r::Array{Float64,3},
                  X_p::Array{Float64,2},
                  X_c::Array{Float64,2})
    (N,N) = size(Y)
    prior(model)+likelihood(model,Y,X_r,X_p,X_c,1:N)
end

function prior(model::ModelState)
    total_prob = sum(normal_logpdf(model.weights))
    total_prob += sum(normal_logpdf(model.beta))
    total_prob += sum(normal_logpdf(model.beta_p))
    total_prob += sum(normal_logpdf(model.beta_r))
    total_prob += sum(normal_logpdf(model.a))
    total_prob += sum(normal_logpdf(model.b))
    total_prob += normal_logpdf(model.c)

    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N = (_2Nm1+1)/2
    total_prob += sum(log(1:N))
    for i = 1:length(tree.nodes)
        total_prob -= log(tree.nodes[i].num_leaves)
        segment_length = (1-gamma)*gamma^(tree.nodes[i].num_ancestors-1)
        total_prob += poisson_logpdf(tree.nodes[i].location, lambda * segment_length)
    end
    total_prob
end

function likelihood(model::ModelState,
                    Y::Array{Int64,2},
                    X_r::Array{Float64,3},
                    X_p::Array{Float64,2},
                    X_c::Array{Float64,2},
                    indices::Array{Int64,1})
    (N,N) = size(Y)
    Z = ConstructZ(tree)
    W = model.weights
    beta = model.beta
    beta_p = model.beta_p
    beta_c = model.beta_r
 
    total_prob = 0 
    for i = indices
        for j = 1:N
            logit_arg = Z[i,:] * W * Z[j,:]' + 
                        beta' * X_r + beta_p' * X_p + beta_c' * X_c + 
                        model.a[i] + model.b[j] + model.c
            if Y[i,j] == 1
                total_prob += -log(1+exp(-logit_arg))
            elseif Y[i,j] == 0
                total_prob += -logit_arg - log(1+exp(-logit_arg))
            end
        end
    end
    total_prob
end

#Local pdfs for sampling

#pdf for updating psi, the tree structure, assumes tree is already pruned at prune_index
function psi_prior_path_logpdf(model::ModelState,
                               path_leaf_index::Int64,
                               pruned_index::Int64)
    (N,N) = size(Y)
    tree = model.tree
    path = GetPath(tree,path_leaf_index)

    psi_probs = prior_path(tree, pruned_index, path_leaf_index)

    mutation_probs = zeros(length(path))
    for i = 1:length(path)
        cur = path[i]

        if i == 1 #leaf node
            poisson_mean_before = model.lambda * gamma^(cur.num_ancestors)
        else
            poisson_mean_before = model.lambda * (1 - gamma) * gamma^(cur.num_ancestors - 1)
        end

        poisson_mean_after = poisson_mean_before*gamma
        mutation_probs[i] = poisson_logpdf(cur.location, poisson_mean_after) -
                            poisson_logpdf(cur.location, poisson_mean_before)
        if i > 1
            mutation_probs[i] += mutation_probs[i-1]
 
            queue = Int64[]

            #perform update for subtree that hasn't been visted yet            
            if cur.children[1].index = path[i-1].index
                enqueue(queue, cur.children[2].index)
            else
                enqueue(queue, cur.children[1].index)
            end

            while length(queue) > 0
                qcur = tree.nodes[pop(queue)]
                if qcur != Nil()
                    enqueue(queue, qcur.children[1].index)
                    enqueue(queue, qcur.children[2].index)
                    poisson_mean_before = model.lambda * (1 - gamma) * gamma^(qcur.num_ancestors - 1)
                    poisson_mean_after = poisson_mean_before*gamma
                    mutation_probs[i] += poisson_logpdf(qcur.location, poisson_mean_after) -
                                         poisson_logpdf(qcur.location, poisson_mean_before)
                end
            end

        end
    end

    max_additional_ancestors = length(path) - 1

    #go through pruned subtree and create hash

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    subtree_probs = zeros(length(path))
    for i = 1:length(path)
        for j = subtree_indices
            cur = tree.nodes[j]
            if j <= N
                poisson_mean = model.lambda * gamma^(cur.num_ancestors + max_additional_ancestors + 1 - i)
            else
                poisson_mean = model.lambda * (1 - gamma) * gamma^(cur.num_ancestors + max_additional_ancestors - i)
            end
            subtree_probs[i] += poisson_logpdf(cur.location, poisson_mean)
        end
    end
    
    psi_probs + mutation_probs + subtree_probs
end

function psi_likelihood_path_logpdf(model::ModelState,
                                    path_leaf_index::Int64,
                                    pruned_index::Int64,
                                    Y::Array{Int64,2},
                                    X_r::Array{Float64,3},
                                    X_p::Array{Float64,2},
                                    X_c::Array{Float64,2})
    (N,N) = size(Y)
    tree = model.tree
    path = GetPath(tree,path_leaf_index)
    


end

# log \prod M_i for all valid grafting points
function prior_tree(tree::Tree,
                    pruned_index::Int)
    N = (length(tree.nodes) + 1) / 2

    # -log \prod M_i for all points below current one, both original and 
    # for grafted trees where the subtree was grafted below the current node
    original_probs = -Inf*ones(2N - 1)
    new_probs = -Inf*ones(2N - 1)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    i = 1
    while contains(subtree_indices, i)
        i += 1
    end
    root = FindRoot(tree, i)

    # stack will contain a leaves to root ordering of nodes in the pruned tree
    stack = GetNodeStack(tree, root.index)

    while length(stack) > 0
        i = pop(stack)
        cur = tree.nodes[i]
        original_probs[i] = -log(max(cur.num_leaves - 1, 1))
        new_probs[i] = -log(cur.num_leaves + tree.nodes[pruned_index].num_leaves - 1)
        if i > N
            for j = 1:2
                child_original_prob = original_probs[cur.children[j].index]
                child_new_prob = new_probs[cur.children[j].index]
                if child_original_prob == -Inf || child_new_prob == -Inf
                    assert(contains(subtree_indices, cur.children[j].index))
                end
                original_probs[i] += child_original_prob
                new_probs[i] += child_new_prob
            end
        end
    end

    new_probs - original_probs
end

# log \prod M_i for a path starting at leaf node_index
function prior_path(tree::Tree,
                    pruned_ind::Int,
                    node_index::Int)
    num_leaves = GetPathNumDescendants(tree, node_index)
    proposed_num_leaves = ProposedNumDescendants(tree,pruned_ind,node_index)
    original_partial_values = cumsum(log(num_leaves-1))
    proposed_partial_values = cumsum(log(proposed_num_leaves-1))
    -(original_partial_values[end]-original_partial_values+proposed_partial_values)
end

#Utility pdfs
function normal_logpdf(x)
    -x.*x/2 - log(sqrt(2pi))
end

function poisson_logpdf(k,lambda)
    -lambda + k*log(lambda) - sum(log(1:k))
end
