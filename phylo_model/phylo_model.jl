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
    ModelState(model.lambda, model.gamma, copy(model.tree), copy(model.Z))
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
function tau(t::TreeNode{T})
    @assert t != Nil()
    p = t.parent
    if p == Nil()
        return p
    end

    c = t 
    while p.children[2] != c
        c = p
        p = p.parent

        if p == Nil()
            return p
        end
    end 

    return p
end
