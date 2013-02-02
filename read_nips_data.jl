using MAT

nips_data = matopen("../data/nips_coauthorship/nips_1-17_matv7.3.mat")

authors_names = read(nips_data["authors_names"]);
aw_counts = read(nips_data["aw_counts"]);
counts = read(nips_data["counts"]);
docs_authors = read(nips_data["docs_authors"]);
docs_names = read(nips_data["docs_names"]);
words = read(nips_data["words"]);

Y = docs_authors' * docs_authors
y_ind_nnz = find(Y)
YY = copy(Y)
YY[y_ind_nnz] = 1
 
authors_numlinks = squeeze(sum(YY, 2))

perm = sortperm(authors_numlinks)
perm = reverse(perm)

top234 = perm[1:234]

Y234 = convert(Array{Int64}, YY[top234, top234])
YY = copy(Y234)

