using Distributions
require("probability_util.jl")

function generate_Z(N, lambda, gam)

    Z = zeros(N,0)

    next_ind = 2N-1
    
    partitions = { ([1:N], 1.0, 2N-1) }

    while length(partitions) > 0
        (points, cur_split, cur_ind) = pop!(partitions)
        if length(points) == 1
            continue
        end
        left_points = Int[]
        right_points = Int[]
        split_point = -1.0
        while length(left_points) == 0 || length(right_points) == 0
            split_point = rand(Uniform(0,1))
            U = rand(Bernoulli(split_point), N)
            left_points = points[find(U[points] .== 0)]
            right_points = points[find(U[points] .== 1)]
            if length(left_points) == 0
                cur_split *= 1-split_point
            elseif length(right_points) == 0
                cur_split *= split_point
            end
        end

        left_ind = next_ind-1
        right_ind = next_ind-2
        next_ind -= 2
       
        left_split = cur_split*split_point
        right_split = cur_split*(1-split_point)

        push!(partitions, (left_points, left_split, left_ind))
        push!(partitions, (right_points, right_split, right_ind))

        tp = cur_split^gam
        tl = length(left_points) == 1 ? 0.0 : left_split^gam
        tr = length(right_points) == 1 ? 0.0 : right_split^gam

        Kl = rand(Poisson(lambda*(tp-tl)))
        Kr = rand(Poisson(lambda*(tp-tr)))

        if Kl > 0    
            Zl = zeros(N, Kl)
            Zl[left_points, :] = 1
            Z = [Z Zl]
        end
        if Kr > 0
            Zr = zeros(N, Kr)
            Zr[right_points, :] = 1
            Z = [Z Zr]
        end

    end

    uof(Z)[1]
end

function generate_synthetic_data(N, lambda, gamma, w_sigma)
    Z = generate_Z(N, lambda, gamma)

    (N,K) = size(Z)
    W = rand(Normal(0,w_sigma),(K,K))

    effects = Z*W*Z'

    pY = exp(broadcast(log_predictive,effects))

    Y = broadcast(p->rand(Bernoulli(p)), pY)

    return Z, W, Y
end
