require("utils/plot_utils.jl")
require("eval/eval_utils.jl")

function eval_phylo_experiments(path, filename_base; p=nothing, offset= 0.0, color="red", label="", kind="diamond")
    filenames = readdir(path)

    delim = contains(filename_base, "pwgs") ? ' ' : ','

    if contains(filename_base, "aldous")
        Ytrue = get_true_clustering()
        depths = [20,30,50,70,100,200,300]
        auprs = zeros(length(depths), 10)

        triu_inds = find(triu(ones(288,288),1))

    end

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

    auprs, p, pts
end

function get_true_clustering()
    Y = zeros(288,288)
    ranges = { 1:108, 109:192, 193:252, 253:288 }
    for k = 1:4
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
