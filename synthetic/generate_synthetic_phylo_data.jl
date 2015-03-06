using Distributions
using HDF5
using JLD

function generate_data(passenger_rate::Float64, read_depth::Int64, mu_r::Float64, mu_v::Float64, kind::ASCIIString)

    N = 3

    if kind == "branch"
        data = Array{Float64}[]

        index = 0
        clusters = [Int64[] for i = 1:N]

        remaining_time = rand(Beta(9,1))

        time = remaining_time 
        n_mutations = rand(Poisson(time * passenger_rate)) + 1
        phi = remaining_time
 
        for n = 1:n_mutations
            index += 1
            d = rand(Binomial(read_depth, 0.8))
            a = rand(Binomial(d, (1-phi)*mu_r + phi*mu_v))
            
            push!(data, [index, index, a, d, mu_r, mu_v])
            push!(clusters[1], index)
        end 

        etas = rand(Dirichlet(4*ones(N-1)))

        for i = 1:N-1
            eta = etas[i]

            time = (1-eta) * remaining_time
            n_mutations = rand(Poisson(time * passenger_rate)) + 1
            
            this_phi = phi*eta

            for n = 1:n_mutations
                index += 1
                d = rand(Binomial(read_depth, 0.8))
                a = rand(Binomial(d, (1-this_phi)*mu_r + this_phi*mu_v))
                
                push!(data, [index, index, a, d, mu_r, mu_v])
                push!(clusters[i+1], index)
            end 
        end

    elseif kind == "chain"
        phi = 1.0

        data = Array{Float64}[]

        index = 0
        clusters = [Int64[] for i = 1:N]
        remaining_time = rand(Beta(9,1)) 
        phi = remaining_time

        for i = 1:N
            eta = rand(Beta(4,4))
            time = remaining_time * (1-eta)
            remaining_time -= time

            n_mutations = rand(Poisson(time * passenger_rate)) + 1
    
            println("phi: $phi")
            for n = 1:n_mutations
                index += 1
                d = rand(Binomial(read_depth, 0.8))
                a = rand(Binomial(d, (1-phi)*mu_r + phi*mu_v))
                
                push!(data, [index, index, a, d, mu_r, mu_v])
                push!(clusters[i], index)
            end 
            phi *= eta

        end
    end

    [data[i][j] for i = 1:index, j = 1:6], clusters

end

function generate_datasets(kind)
    for i = 1:8
        A, clusters = generate_data(1000.0, 50, 0.999, 0.5, kind)
        AA = Array(Any, size(A,1)+1, size(A,2))

        AA[2:end,:] = A
        AA[1,:] = ["id", "gene", "a", "d", "mu_r", "mu_v"]

        f = open("betasplit_$(kind)_$i.csv", "w"); writedlm(f,AA,"\t"); close(f)
        save("betasplit_$(kind)_clusters_$i.jld", "clusters", clusters)
    end
end

