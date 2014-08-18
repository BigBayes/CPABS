
function save(filename, array) # can also be a array
    stream = open(filename, "w")

    array_size::Int64 = length(array)

    write(stream,array_size)
    nd = zeros(Int64, array_size)
    for i = 1:array_size
        nd[i] = ndims(array[i])
    end
    write(stream, nd)


    lengths = zeros(Int64, sum(nd))
    count = 0
    for i = 1:array_size
        for j = 1:nd[i]
            count += 1
            lengths[count] = size(array[i])[j]
        end
    end
    write(stream, lengths)

    for i = 1:array_size
        a = convert(Array{Float64}, array[i])
        write(stream, a)    
    end
 
    close(stream)
end

function restore(filename)
    stream = open(filename, "r")
    array_size = read(stream,Int64)
    nds = read(stream, Int64, array_size)
    lengths = read(stream, Int64, sum(nds))

    sizes = [zeros(Int64,0) for i = 1:array_size]
    count = 0
    for i = 1:array_size
        sz = sizes[i]
        for j = 1:nds[i]
            count += 1
            push!(sz,lengths[count]) 
        end
    end

    arrays = [ read(stream, Float64, sizes[i]...) for i = 1:array_size]

    close(stream)
    arrays
end

function restore_in_range(filename_base,
                          index_range)

    [restore("$filename_base$i.jla") for i = index_range]
end



