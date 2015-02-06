require("phylo_model/phylo_model.jl")

function read_phylosub_data(filename)
    if contains(filename, "emptysims")
        data = readdlm("../data/phylosub/$filename", '\t', header=true)
    elseif contains(filename, "aldous")
        data = readdlm("../data/phylosub/$filename", '\t', header=true)
    elseif contains(filename, "phylospan")
        data = readdlm("../data/phylosub/$filename", '\t', header=true)
    elseif contains(filename, "betasplit")
        data = readdlm("../data/phylosub/$filename", '\t', header=true)
    else
        data = readdlm("../data/phylosub/$filename", ';', header=true)
    end
    header = data[2]
    name_index = find(header .== "gene")[1]
    a_index = find(header .== "a")[1]
    d_index = find(header .== "d")[1]
    mu_r_index = find(header .== "mu_r")[1]
    mu_v_index = find(header .== "mu_v")[1]

    num_entries, _7 = size(data[1])

    if typeof(data[1][1,a_index]) == Float64
        S = 1
    else
        a = readdlm(IOBuffer(convert(ASCIIString, data[1][1,a_index])), ',')
        S = length(a)
    end

    names = Array(ASCIIString, num_entries)
    AA = zeros(num_entries, S)
    DD = zeros(num_entries, S)
    mu_r = zeros(num_entries)
    mu_v = zeros(num_entries) 

    for i = 1:num_entries
        names[i] = lowercase(string(data[1][i,name_index]))
        if typeof(data[1][i,a_index]) == Float64
            a = data[1][i,a_index]
            d = data[1][i,d_index]
            AA[i,:] = a
            DD[i,:] = d
        else
            a = readdlm(IOBuffer(convert(ASCIIString, data[1][i,a_index])), ',')
            d = readdlm(IOBuffer(convert(ASCIIString, data[1][i,d_index])), ',')
            AA[i,:] = a[:]
            DD[i,:] = d[:]
        end

        mu_r[i] = data[1][i,mu_r_index]
        mu_v[i] = data[1][i,mu_v_index]
    end
    return (AA, DD, mu_r, mu_v, names)
end

function read_multilocus_data(filename, name2index::Dict{ASCIIString, Int64})
    data, headers = readdlm("../data/phylosub/$filename", '\t', header=true)

    num_entries, _ = size(data)

    paired_reads = zeros(num_entries, 9)
    name_pairs = Array(ASCIIString, (num_entries, 2))
    
    for i = 1:num_entries
        paired_reads[i,1:2] = [name2index[data[i,1]], name2index[data[i,2]]]
        paired_reads[i,3:end] = float64(data[i,3:end])
        paired_reads[i,4] += 1
        #paired_reads[i,5] = 0.5
    end

    paired_reads
end

function constructDataState(filename; multilocus_filename=nothing)
    (AA, DD, mu_r, mu_v, names) = read_phylosub_data(filename)
    paired_reads = zeros(0,9)
    if multilocus_filename != nothing
        name2index = Dict{ASCIIString, Int64}()

        for i = 1:length(names)
            name2index[names[i]] = i
        end

        paired_reads = read_multilocus_data(multilocus_filename, name2index)
    end

    DataState(AA, DD, mu_r, mu_v, paired_reads, names)
end
