require("utils/plot_utils.jl")
require("eval/eval_utils.jl")

function eval_phylo_experiments(path, filename_base; p=nothing, offset= 0.0, color="red")
    filenames = readdir(path)

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

            Ypred = readdlm("$path/$fname", ',')

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
    add(p, Points(depths+offset, meds, color=color, size=0.5))
    display(p)

    auprs, p
end

function get_true_clustering()
    Y = zeros(288,288)
    ranges = { 1:108, 109:192, 193:252, 253:288 }
    for k = 1:4
        Y[ranges[k], ranges[k]] = 1
    end
    Y
end
