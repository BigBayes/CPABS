
function read_phylosub_data(filename)
    if contains(filename, "emptysims")
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
        names[i] = data[1][i,name_index]
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
