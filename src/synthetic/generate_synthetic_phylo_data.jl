#using Distributions
#using HDF5
#using JLD

function generate_data(N::Int64, passenger_rate::Float64, read_depth::Int64, mu_r::Float64, mu_v::Float64, kind::ASCIIString; 
                       num_mutations = nothing, fixed_depth::Bool = false)

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
        phis = zeros(N)
        for i = 1:N
            phis[i] = phi
            eta = rand(Beta(4,4))
            time = remaining_time * (1-eta)
            remaining_time -= time

            n_mutations = num_mutations == nothing ? rand(Poisson(time * passenger_rate)) + 1 : num_mutations
    
            println("phi: $phi")
            for n = 1:n_mutations
                index += 1
                d = fixed_depth ? read_depth : rand(Binomial(read_depth, 0.8))
                a = rand(Binomial(d, (1-phi)*mu_r + phi*mu_v))
                
                push!(data, [index, index, a, d, mu_r, mu_v])
                push!(clusters[i], index)
            end 
            phi *= eta

        end
    end

    [data[i][j] for i = 1:index, j = 1:6], clusters, phis

end

function generate_datasets(kind)
    for i = 1:40
        A, clusters, phis = generate_data(3, 200.0, 50, 0.999, 0.5, kind)
        AA = Array(Any, size(A,1)+1, size(A,2))

        AA[2:end,:] = A
        AA[1,:] = ["id", "gene", "a", "d", "mu_r", "mu_v"]

        f = open("betasplit_$(kind)_$i.csv", "w"); writedlm(f,AA,"\t"); close(f)
        save("betasplit_$(kind)_clusters_$i.jld", "clusters", clusters)
    end
end

function generate_betasplit_phylo_datasets()
    for N = 3:5
        for depth = [50,70,100]
            for num_mutations = [10,25,100]
                for i = 1:8
                    A, clusters, phi = generate_data(N, 1.0, depth, 0.999, 0.5, "chain", 
                                                     num_mutations=num_mutations, fixed_depth=true)
                    AA = Array(Any, size(A,1)+1, size(A,2))

                    AA[2:end,:] = A
                    AA[1,:] = ["id", "gene", "a", "d", "mu_r", "mu_v"]

                    f = open("betasplit_phylo_chain_$(N)_$(depth)_$(num_mutations)_$i.csv", "w"); writedlm(f,AA,"\t"); close(f)
                    @save "gt_betasplit_phylo_chain_$(N)_$(depth)_$(num_mutations)_$i.jld" clusters phi
                end
            end
        end

    end
end
