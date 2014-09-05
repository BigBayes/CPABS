require("tree.jl")

import Base.copy
import Base.ref
import Base.assign

type ModelState
    lambda::Float64
    gamma::Float64
    alpha::Float64

    tree::Tree{Vector{Float64}} # Tree state holds the \eta variables (one eta variable for each observed sample)
    Z::Vector{Int64} # Assignment of datapoints to nodes
end

function copy(model::ModelState)
    ModelState(model.lambda, model.gamma, model.alpha, copy(model.tree), copy(model.Z))
end

type ModelSpecification
    rrj_jump_probabilities::Array{Float64} #assumes L \in {k-1,k,k+1}

    debug::Bool
    verbose::Bool
    plot::Bool
end

copy(ms::ModelSpecification) = ModelSpecification(ms.rrj_jump_probabilities, ms.debug,
                                                  ms.verbose, ms.plot)

type DataState
    reference_counts::Matrix{Float64}
    total_counts::Matrix{Float64}
    mu_r::Vector{Float64}
    mu_v::Vector{Float64}
end

# Find ancestor to whose population t belongs
# As we assume the right children are the "new" subpopulations,
# tau will be the most recent ancestor such that the path from tau to
# t contains tau's right child
function tau(t::TreeNode)
    @assert t != Nil()
    p = t.parent
    if p == Nil()
        return p
    end

    c = t 
    while p.children[1] != c
        c = p
        p = p.parent

        if p == Nil()
            return p
        end
    end 

    return p
end

# Like tau, but returns the full path from t to tau(t)
function tau_path(t::TreeNode)
    @assert t != Nil()
    path = Array(TreeNode,0)
    push!(path,t)
    p = t.parent
    if p == Nil()
        push!(path,p)
        return path
    end

    c = t 
    while p.children[1] != c
        c = p
        p = p.parent
        push!(path,p)

        if p == Nil()
            return path
        end
    end 

    return path

end

function compute_times(model::ModelState)
    tree = model.tree
    gam = model.gamma
    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    t = ones(2N-1)

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            t[i] = t[parent.index]*(cur.rhot*cur_mu_prop)^gam

        else
            t[i] = cur.rhot^gam
        end

    end

    return t
end

function compute_taus(model::ModelState)
    tree = model.tree

    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    Tau = zeros(2N-1)
    indices = GetLeafToRootOrdering(tree, root.index)

    for i = reverse(indices)
        cur = tree.nodes[i]
        if i == root.index
            Tau[i] = 0
        end

        if i > N
            l = cur.children[2].index
            r = cur.children[1].index

            Tau[l] = Tau[i]
            Tau[r] = i 
        end 
    end

    return Tau
end

function compute_phis(model::ModelState)
    tree = model.tree

    S = length(tree.nodes[end].state)
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    Tau = zeros(2N-1)
    indices = GetLeafToRootOrdering(tree, root.index)

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)[1]

            for s = 1:S
                # the eta variable held by a node is the eta for the right child
                eta_self = self_direction == 1 ? parent.state[s] : (1 - parent.state[s])
                B[i,s] = B[parent.index]*eta_self
            end
        else
            B[i,:] = 1.0
        end

    end

    return B
end

function model2array(model::ModelState)
    tree2array(model.tree, model.gamma)
end
