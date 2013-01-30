
function save(filename, array) # can also be a array
    stream = open(filename, "w")

    array_size::Int64 = length(array)
    write(stream, array_size)

    lengths = [length(array[i]) for i = 1:array_size]
    write(stream, lengths)

    for i = 1:array_size
        a = convert(Array{Float64,1}, array[i])
        write(stream, a)    
    end
 
    close(stream)
end

function restore(filename)
    stream = open(filename, "r")
    array_size = read(stream,Int64)
    lengths = read(stream, Int64, array_size)

    arrays = [ read(stream, Float64, lengths[i]) for i = 1:array_size]

    close(stream)
    arrays
end
