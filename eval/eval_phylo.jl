require("utils/plot_utils.jl")
require("eval/eval_utils.jl")

function eval_phylo_experiments(path, filename_base; p=nothing, offset= 0.0, color="red", label="", kind="diamond", alpha=1)
    filenames = readdir(path)

    delim = contains(filename_base, "pwgs") ? ' ' : ','

    if contains(filename_base, "aldous")
        Ytrue = get_true_clustering()
        depths = [20,30,50,70,100,200,300]
        auprs = zeros(length(depths), 10)

        triu_inds = find(triu(ones(288,288),1))


        for fname in filenames
            if contains(fname, filename_base)
                m = match(r"\.([0-9]+)\.([0-9]+)\.csv", fname)
                D = int(m.captures[1])
                count = int(m.captures[2])

                D_index = find(D .== depths)[1]

                Ypred = readdlm("$path/$fname", delim)
                auprs[D_index, count] = aupr(Ypred[triu_inds], Ytrue[triu_inds])
            end
        end
        meds = median(auprs,2)
        mins = minimum(auprs,2)
        maxs = maximum(auprs,2)

        if p == nothing
            p = FramedPlot()
        end

        ebars = Winston.ErrorBarsY(depths+offset, mins, maxs, color=color)
        add(p, ebars)
        pts = Points(depths+offset, meds, color=color, size=0.5, kind=kind)
        setattr(pts, label=label)
        add(p, pts)
        display(p)

        return auprs, p, pts

    elseif contains(filename_base, "emptysim")

        n_clusters = [2,3,4,5]
        depths = [20,30,50,70,100,200,300]
        n_mutations = [5,10,25,50,100,200,500]


        auprs = zeros(length(n_clusters), length(depths), length(n_mutations))

        for fname in filenames
            if contains(fname, filename_base)
                m = match(r"\.1\.0\.([0-9]+)\.([0-9]+)\.([0-9]+)\.csv", fname)
                C = int(m.captures[1])
                D = int(m.captures[2])
                N = int(m.captures[3])
       
                if N > 500
                    continue
                end

                triu_inds = find(triu(ones(C*N,C*N),1))

                C_index = find(C .== n_clusters)[1]
                D_index = find(D .== depths)[1]
                N_index = find(N .== n_mutations)[1]

                Ypred = readdlm("$path/$fname", delim)
                Ytrue = get_true_clustering_emptysims(C, N)

                auprs[C_index, D_index, N_index] = aupr(Ypred[triu_inds], Ytrue[triu_inds])
            end
        end

        return auprs
    end
end

function get_true_clustering()
    Y = zeros(288,288)
    ranges = { 1:108, 109:192, 193:252, 253:288 }
    for k = 1:4
        Y[ranges[k], ranges[k]] = 1
    end
    Y
end

function get_true_clustering_emptysims(C, N)
    base_range = 1:N
    ranges = { base_range + i*N for i = 0:C-1}

    Y = zeros(C*N,C*N)
    for k = 1:C
        Y[ranges[k], ranges[k]] = 1
    end
    Y
end

function plot_paired_scatter(auprs1, auprs2, label1, label2, row_index, title)

    color1 = "red"
    color2 = "blue"

    kind1 = "filled circle"
    kind2 = "filled square"

    p = FramedPlot()

    c1 = Points([1:10], auprs1[row_index, :], kind=kind1, color=color1)
    c2 = Points([1:10], auprs2[row_index, :], kind=kind2, color=color2)
   
    setattr(c1, label=label1)
    setattr(c2, label=label2)

    add(p,c1)
    add(p,c2)

    min_aupr = minimum([auprs1[row_index,:] , auprs2[row_index, :]])

    max_aupr = maximum([auprs1[row_index,:] , auprs2[row_index, :]])
    
    setattr(p, yrange=( min_aupr - 0.15, max_aupr+0.03))

    l = Legend(0.65, 0.2, {c1, c2})

    add(p, l)

    setattr(p, xlabel="Trial")
    setattr(p, ylabel="AUPR")
    setattr(p, title=title)
    display(p)

    p    
end

function read_phylosub_results()
    A, header = readdlm("../results/emsim.results.txt",'\t', header=true)

    n_clusters = [2,3,4,5]
    depths = [20,30,50,70,100,200,300]
    n_mutations = [5,10,25,50,100,200,500]

    auprs = zeros(length(n_clusters), length(depths), length(n_mutations))

    for i = 1:size(A,1)
        C = int(A[i,1])-1
        D = int(A[i,2])
        N = int(A[i,3])

        if N > 500
            continue
        end

        C_index = find(C .== n_clusters)[1]
        D_index = find(D .== depths)[1]
        N_index = find(N .== n_mutations)[1]

        auprs[C_index, D_index, N_index] = float(A[i,4])
    end

    auprs
end
