
type WangLandauState
    energy_boundaries::Vector{Float64} # e.g. [-1000, -10, Inf]
    K_boundaries::Vector{Float64} # model indices that define the boundaries for the partition of the models, 
                                # eg [3,5,6,7,Inf] defines the partition 1-2, 3-4, 5, 6, 7-Inf
    partition_function::Matrix{Float64} # (length(K_boundaries) x (length(energy_boundaries)+1)
    histogram::Matrix{Int64} # (max_K+1) x (length(energy_boundaries)+1)

    initial_stepsize::Float64
    stepsize_epoch::Int64 
    stepsize::Function
    flat_histogram_test::Function
end

function WangLandauState(energy_boundaries::Vector{Float64},
                         K_boundaries::Vector{Float64}, 
                         initial_stepsize::Float64,
                         histogram_test_frequency_ratio::Float64)

    if energy_boundaries[end] != Inf
        push!(energy_boundaries, Inf)
    end

    if K_boundaries[end] != Inf
        push!(K_boundaries, Inf)
    end

    partition_function = zeros(length(K_boundaries), length(energy_boundaries))
    histogram = zeros(Int64, length(K_boundaries), length(energy_boundaries))

    stepsize_epoch = 1
    stepsize_function = wang_landau_stepsize
    flat_histogram_test = gen_wang_landau_histogram_test(histogram_test_frequency_ratio)

    WangLandauState(energy_boundaries, K_boundaries, partition_function, histogram, initial_stepsize,
                    stepsize_epoch, stepsize_function, flat_histogram_test)
end

# WL without any partitions
function WangLandauState()
    WangLandauState([Inf], [Inf], 0.0, 0.0) 
end    

function copy(wl_state::WangLandauState)
    WangLandauState(copy(wl_state.energy_boundaries), copy(wl_state.K_boundaries), 
                    copy(wl_state.partition_function), copy(wl_state.histogram),
                    wl_state.initial_stepsize, wl_state.stepsize_epoch, 
                    wl_state.stepsize, wl_state.flat_histogram_test)
end

function stepsize(wl_state::WangLandauState)
    wl_state.stepsize(wl_state.initial_stepsize, wl_state.stepsize_epoch)
end

function wang_landau_stepsize(f0::Float64, epoch::Int64)
    f0/(2*epoch)
end

function gen_wang_landau_histogram_test(frequency_ratio::Float64)

    function wang_landau_histogram_test(histogram::Matrix{Int64})
        hmean = mean(histogram)
        hmin = minimum(histogram)

        hmin > frequency_ratio * hmean
    end

    wang_landau_histogram_test
end

function update_partition_function(wl_state::WangLandauState, K::Int64, energy::Float64)

    mod_factor = stepsize(wl_state)
    energy_boundaries = wl_state.energy_boundaries
    K_boundaries = wl_state.K_boundaries

    K_bin = find(K .< K_boundaries)[1] 
    energy_bin = find(energy .< energy_boundaries)[1]

    wl_state.partition_function[K_bin, energy_bin] += mod_factor
    wl_state.histogram[K_bin, energy_bin] += 1

    if wl_state.flat_histogram_test(wl_state.histogram)
        wl_state.stepsize_epoch += 1
        wl_state.histogram = zeros(size(wl_state.histogram))
    end
 
end

function get_partition_function(wl_state::WangLandauState, K::Int64, energy::Float64)
    @assert energy < Inf
    energy_boundaries = wl_state.energy_boundaries
    K_boundaries = wl_state.K_boundaries

    K_bin = find(K .< K_boundaries)[1] 
    energy_bin = find(energy .< energy_boundaries)[1]
    wl_state.partition_function[K_bin, energy_bin]
end

function get_bin(wl_state::WangLandauState, K::Int64, energy::Float64)
    @assert energy < Inf
    energy_boundaries = wl_state.energy_boundaries
    K_boundaries = wl_state.K_boundaries

    K_bin = find(K .< K_boundaries)[1] 
    energy_bin = find(energy .< energy_boundaries)[1]

    K_bin, energy_bin
end
