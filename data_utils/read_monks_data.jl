
Y1 = readcsv("../data/monks/sampson_agent_agent[SAMPLK1].csv")
Y2 = readcsv("../data/monks/sampson_agent_agent[SAMPLK2].csv")
Y3 = readcsv("../data/monks/sampson_agent_agent[SAMPLK3].csv")

y_ind_nnz = find(Y1)
YY1 = convert(Array{Float64},Y1)
YY1[y_ind_nnz] = 1

y_ind_nnz = find(Y2)
YY2 = convert(Array{Float64},Y2)
YY2[y_ind_nnz] = 1

y_ind_nnz = find(Y3)
YY3 = convert(Array{Float64},Y3)
YY3[y_ind_nnz] = 1

YY = Array(Array{Float64,2},3)
YY[1] = YY1
YY[2] = YY2
YY[3] = YY3
