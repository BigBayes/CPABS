
type WangLandauState
    energy_boundaries::Vector{Float64} # e.g. [-1000, -10, Inf]
    partition_function::Matrix{Float64} # (max_K+1) x (length(energy_boundaries)+1)
    histogram::Matrix{Int64} # (max_K+1) x (length(energy_boundaries)+1)
   
    max_K::Int64 # maximum model index to have its own partition_function row
    initial_stepsize::Float64
    stepsize_epoch::Int64 
    stepsize::Function
    flat_histogram_test::Function
end

function WangLandauState(energy_boundaries::Vector{Float64},
                         max_K::Int64, 
                         initial_stepsize::Float64,
                         histogram_test_frequency_ratio::Float64)

    if energy_boundaries[end] != Inf
        push!(energy_boundaries, Inf)
    end

    partition_function = zeros(max_K+1, length(energy_boundaries))
    histogram = zeros(Int64, max_K+1, length(energy_boundaries))

    stepsize_epoch = 1
    stepsize_function = wang_landau_stepsize
    flat_histogram_test = gen_wang_landau_histogram_test(histogram_test_frequency_ratio)

    WangLandauState(energy_boundaries, partition_function, histogram, max_K, initial_stepsize,
                    stepsize_epoch, stepsize_function, flat_histogram_test)
end
    

function copy(wl_state::WangLandauState)
    WangLandauState(copy(wl_state.energy_boundaries), copy(wl_state.partition_function), copy(wl_state.histogram),
                    wl_state.modification_factor, wl_state.max_K, wl_state.stepsize_epoch, wl_state.initial_stepsize,
                    wl_state.stepsize_update, wl_state.flat_histogram_test)
end

function stepsize(wl_state::WangLandauState)
    wl_state.stepsize(wl_state.initial_stepsize, wl_state.stepsize_epoch)
end

function wang_landau_stepsize(f0::Float64, epoch::Int64)
    log(1+f0)/epoch
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

    K = K > max_K ? max_K+1 : K
    energy_bin = find(energy .< energy_boundaries)[1]

    wl_state.partition_function[K, energy_bin] += mod_factor
    wl_state.histogram[K, energy_bin] += 1

    if wl_state.flat_histogram_test(wl_state.histogram)
        stepsize_epoch += 1
        wl_state.histogram = zeros(size(wl_state.histogram))
    end
 
end

function get_partition_function(wl_state::WangLandauState, K::Int64, energy::Float64)
    K = K > max_K ? max_K+1 : K
    energy_bin = find(energy .< energy_boundaries)[1]
    wl_state.partition_function[K, energy_bin]
end

