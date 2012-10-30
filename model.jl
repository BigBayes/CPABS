type ModelState
    lambda::Float64
    gamma::Float64
    tree::Tree{Int64}
    weights::Array{Float64,2}
    beta::Array{Float64,1}
    beta_p::Array{Float64,1}
    beta_c::Array{Float64,1}
    a::Array{Float64,1}
    b::Array{Float64,1}
    c::Float64
end
