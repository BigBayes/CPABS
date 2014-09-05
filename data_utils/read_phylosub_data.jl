if !isdefined(:filename)
    filename = "CLL077.csv"
end

data = readdlm("../data/phylosub/$filename", ';', header=true)

num_entries, _7 = size(data[1])

a = readdlm(IOBuffer(convert(ASCIIString, data[1][1,2])), ',')
S = length(a)

names = Array(ASCIIString, num_entries)
AA = zeros(num_entries, S)
DD = zeros(num_entries, S)
mu_r = zeros(num_entries)
mu_v = zeros(num_entries) 

for i = 1:num_entries
    names[i] = data[1][i,1]
    a = readdlm(IOBuffer(convert(ASCIIString, data[1][i,2])), ',')
    d = readdlm(IOBuffer(convert(ASCIIString, data[1][i,3])), ',')
    AA[i,:] = a[:]
    DD[i,:] = d[:]

    mu_r[i] = data[1][i,4]
    mu_v[i] = data[1][i,6]
end
