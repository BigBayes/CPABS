
function auc(logprobs::Vector{Float64}, targets::Vector{Float64})
    N = length(targets)
    Npos = sum(targets)
    Nneg = N - Npos

    s_ind = sortperm(logprobs)

    t_ind = [0, find(targets[s_ind] .== 1)]
    t_diff = t_ind[2:end] .- t_ind[1:end-1] .- 1 # number of 0s hit before seeing a 1
    sum(t_diff .* [Npos:-1:1]) / (Npos*Nneg)
end

function aupr(predictions::Vector{Float64}, targets::Vector{Float64})
    N = length(targets)
    Npos = sum(targets)
    Nneg = N - Npos
   
    s_ind = sortperm(predictions)

    s_targets = targets[s_ind]
    reverse!(s_targets)

    precisions = cumsum(s_targets) ./ [1:N]
    recalls = cumsum(s_targets) ./ Npos

    s_ind = sortperm(recalls)
    precisions = [1.0, precisions[s_ind]]
    recalls = [0.0, recalls[s_ind]]

    # trapezoidal area
    sum(0.5*(precisions[2:end] .+ precisions[1:end-1]) .* (recalls[2:end] .- recalls[1:end-1]) )
end

function get_vectors(chain, variable_name::ASCIIString, vectorize::Function)

    v1 = vectorize(chain[1][1][variable_name])
    V = zeros(length(chain), length(v1))

    for i = 1:length(chain)
        x = chain[i][1][variable_name]
        V[i,:] = vectorize(x)
    end

    V
end
