
Y = readcsv("../data/monks/sampson_agent_agent[SAMPLK3].csv")

y_ind_nnz = find(Y)
YY = convert(Array{Int64},Y)
YY[y_ind_nnz] = 1
 
