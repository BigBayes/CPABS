require("tree.jl")

import Base.copy
import Base.ref
import Base.assign

type ModelState
    lambda::Float64
    gamma::Float64
    tree::Tree{Float64} # Tree state holds the \eta variables, represented as gamma draws
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
    read_frequencies::Vector{Float64}
end
