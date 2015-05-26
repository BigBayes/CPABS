
function auc(logprobs::Vector{Float64}, targets::Vector{Float64})
    N = length(targets)
    Npos = sum(targets)
    Nneg = N - Npos

    s_ind = sortperm(logprobs)

    t_ind = [0, find(targets[s_ind] .== 1)]
    t_diff = t_ind[2:end] .- t_ind[1:end-1] .- 1 # number of 0s hit before seeing a 1
    sum(t_diff .* [Npos:-1:1]) / (Npos*Nneg)
end

function aupr(predictions::Vector{Float64}, targets::Vector{Float64}; return_curves=false)
    N = length(targets)
    Npos = sum(targets)

    if Npos == 0
        return 0.0
    end

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
    A = sum(0.5*(precisions[2:end] .+ precisions[1:end-1]) .* (recalls[2:end] .- recalls[1:end-1]) )
    if return_curves
        return (A, precisions, recalls)
    else
        return A
    end
end

function adjusted_mutual_information(predicted_clusters, true_clusters)
    K1 = length(predicted_clusters)
    K2 = length(true_clusters)


    n1 = [length(A) for A in predicted_clusters]
    n2 = [length(A) for A in true_clusters]

    @assert all(n1 .> 0)
    @assert all(n2 .> 0)

    n = [length(intersect(U,V)) for U in predicted_clusters, V in true_clusters]
    N = sum(n1) 
    
    EMI = 0.0
    MI = 0.0

    for i = 1:K1
        for j = 1:K2
            ai = n1[i]
            bj = n2[j]
            for nij = max(1, ai + bj - N):min(ai, bj)
                fterm = lfact(ai) + lfact(bj) + lfact(N-ai) + lfact(N-bj) - 
                        lfact(N) - lfact(nij) - lfact(ai-nij) - lfact(bj-nij) - lfact(N-ai-bj+nij)
                EMI += nij/N * log(N*nij/(ai*bj)) * exp(fterm)
            end

            p12 = n[i,j]/N
            p1 = n1[i]/N
            p2 = n2[j]/N

            MI += p12 > 0 ? p12*log(p12/(p1*p2)) : 0.0
            
        end
    end 

    H1 = - sum([ ni/N * log(ni/N) for ni in n1])
    H2 = - sum([ ni/N * log(ni/N) for ni in n2])

    return (MI - EMI)/(max(H1,H2) - EMI)
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
