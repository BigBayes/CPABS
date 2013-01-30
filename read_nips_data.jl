load("matio.jl")

nips_data = MatIO.matopen("../data/nips_coauthorship/nips_1-17_matv7.3.mat")

authors_names = read(nips_data["authors_names"]);
aw_counts = read(nips_data["aw_counts"]);
counts = read(nips_data["counts"]);
docs_authors = read(nips_data["docs_authors"]);
docs_names = read(nips_data["docs_names"]);
words = read(nips_data["words"]);

Y = docs_authors' * docs_authors
authors_numlinks = squeeze(sum(Y, 2))

sorted,perm = sortperm(authors_numlinks)
perm = reverse(perm)

top234 = perm[1:234]

Y234 = convert(Array{Int64}, Y[top234, top234])
y_ind_nnz = find(Y234) 
YY = copy(Y234)
YY[y_ind_nnz] = 1


