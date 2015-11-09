using Winston
import Winston.Interval
import Winston.imagesc

using Compose
using Color
#function imagegs{T<:Real}(xrange::Interval, yrange::Interval, data::AbstractArray{T,2}, clims::Interval; plot=true)
#    p = FramedPlot()
#    setattr(p, "xrange", xrange)
#    setattr(p, "yrange", reverse(yrange))
#    if clims[1] == clims[2] # avoid divide-by-0
#        clims = (clims[1], clims[1]+1)
#    end
#
#    img = Winston.data2rgb(data, clims, Winston.GrayColormap())
#    add(p, Image(xrange, reverse(yrange), img))
#
#    if plot
#        Winston.tk(p)
#    else
#        return p
#    end
#end
#
#imagesc(xrange, yrange, data) = imagesc(xrange, yrange, data, (minimum(data), minimum(data) == maximum(data) ? maximum(data)+1 : maximum(data)))
#imagegs(data; plot=true) = ((h, w) = size(data); imagegs((0,w), (0,h), max(data)-data, plot=plot))

# Upper order form -- push zeros upwards, giving leftwards columns priority 
function uof(Z::Matrix)
    (N,K) = size(Z)

    function lt(v1, v2)
        for i = 1:length(v1)
            if v1[i] < v2[i]
                return true
            elseif v2[i] < v1[i]
                return false
            end
        end
        return false
    end

    Zarray = Array(Vector{Int},N)
    for i = 1:N
        Zarray[i] = Z[i,:][:]
    end

    perm = sortperm(Zarray, lt=lt)

    Znew = Z[perm,:]   
    Znew, perm
end

function dendrogram(Z, U; plot=true, labels=nothing, leaf_times=nothing, sorted_inds=nothing, annotations=nothing)
    U = convert(Array{Int64}, U)
    Nm1, _ = size(Z)
    N = Nm1+1
    # Sort in depth-first-search order so the dendrogram doesn't cross itself 
    times = zeros(2N-1)
    if leaf_times!=nothing
        times[1:N] = leaf_times
    end
    locations = zeros(2N-1)
    parents = zeros(2N-1)

    visited_stack = [2N-1] 

    leaf_order = Int[]
    while length(visited_stack) > 0
        ind = pop!(visited_stack)

        if ind <= N
            push!(leaf_order, ind)
            continue
        end

        l = Z[ind-N,1]
        r = Z[ind-N,2]
        times[ind] = Z[ind-N,3]

        parents[l] = ind
        parents[r] = ind
        
        push!(visited_stack, r) 
        push!(visited_stack, l)
    end

    locations[leaf_order] = [0:N-1]/(N-1)

    # Traverse in breadth first order so get leaves-to-root ordering of internal nodes
    stack = Int[]
    queue = [2N-1]
    while length(queue) > 0
        ind = pop!(queue)
        l = Z[ind-N,1]
        r = Z[ind-N,2]
       
        if l > N
            unshift!(queue, l)
        end
        if r > N
            unshift!(queue, r)
        end

        push!(stack, ind) 
    end
    print_order = reverse(stack)
  
    mutations_x = Float64[]
    mutations_y = Float64[]
 
    p = FramedPlot()
    for n = print_order
        l = Z[n-N,1]
        r = Z[n-N,2]
        l_loc = locations[l]
        r_loc = locations[r]
        locations[n] = (l_loc + r_loc)/2

        l_t = times[l]
        r_t = times[r]
        n_t = times[n]

        x_points = [l_loc, l_loc, r_loc, r_loc]
        y_points = [l_t, n_t, n_t, r_t]

        add(p, Curve(x_points, y_points, "color", "blue"))

        if sorted_inds != nothing
            original_l = sorted_inds[int(l)]
            original_r = sorted_inds[int(r)]
            add(p, Winston.DataLabel(l_loc, l_t-0.02, "$original_l", color="red", size=0.3))
            add(p, Winston.DataLabel(r_loc, r_t-0.02, "$original_r", color="red", size=0.3))

            if annotations != nothing
                add(p, Winston.DataLabel(l_loc+0.06, l_t+0.05, annotations[original_l], color="black", size=0.3))
                add(p, Winston.DataLabel(r_loc+0.06, r_t+0.05, annotations[original_r], color="black", size=0.3))
            end
        end

        if U[l] > 0
            l_x, l_y = get_mutation_plot_locations(U[l], l_loc, l_t, n_t)
            append!(mutations_x, l_x)
            append!(mutations_y, l_y)
        end
        if U[r] > 0
            r_x, r_y = get_mutation_plot_locations(U[r], r_loc, r_t, n_t)
            append!(mutations_x, r_x)
            append!(mutations_y, r_y)
        end

        if n == 2N-1
            if U[n] > 0
                original_n = sorted_inds[int(n)]
                n_x, n_y = get_mutation_plot_locations(U[n], locations[n], n_t, 1.0)
                append!(mutations_x, n_x)
                append!(mutations_y, n_y)
                add(p, Winston.DataLabel((l_loc+r_loc)/2, n_t-0.02, "$original_n", color="red", size=0.3))
                if annotations != nothing
                    add(p, Winston.DataLabel((l_loc+r_loc)/2+0.06, n_t+0.05, annotations[original_n], color="black", size=0.3))
                end
            end
        end
    end

    add(p, Points(mutations_x, mutations_y, kind="diamond", size=0.3))

    if labels != nothing
        setattr(p.x1, "ticks", locations[1:N])
        setattr(p.x1, "ticklabels", labels)
    end

    if plot
        display(p)
    end
    p
end

# T is the adjacency matrix of the tree to plot
# Z is a dictionary containing a string to display for each cluster
# Heavily borrowed from Compose.introspect()
function plot_subclonal_hierarchy(T, Z; labels=nothing)
    # Find root
    root = 0
    for c = 1:size(T,1)
        if sum(T[:,c]) == 0 && sum(T[c,:]) > 0
            @assert root == 0
            root = c
        end
    end

    positions = Dict{Int64, (Float64, Float64)}()
    queue = Any[]
    push!(queue, (root, 0.0))

    figure_width = Z == nothing ? 10cm : 20cm

    Compose.set_default_graphic_size(figure_width, 10cm)
    figsize = 10mm
    figs = compose!(context(), stroke("#333"), linewidth(0.5mm))
    #p = FramedPlot()
    #add(p, Winston.DataLabel( 0.0, 0.0, label, color="black", size=2.5))
    label = labels == nothing ? "$root" : labels[root]

    cur_x = -1.0
    depth = -Inf
    max_x = -Inf
    while !isempty(queue)

        cur, cur_y = shift!(queue)

        if cur_y > depth
            cur_x = -1.0
        end
        cur_x += 1
        children = find(T[cur,:])
        L = length(children)

        label = labels == nothing ? "$cur" : labels[cur]

        fig = context(cur_x, cur_y)

        compose!(fig, circle(0.5, 0.5, figsize/2), fill("#333"))
        label_fig = context(cur_x, cur_y, order=1)
        compose!(label_fig, text(0.5, 0.5, label, hcenter, vcenter, Rotation()), stroke("white"), fill("white"), fontsize(8))

        positions[cur] =  (cur_x+0.5, cur_y+0.5)

        for c in children
            push!(queue, (c,cur_y+1))
        end

        max_x = max(max_x, cur_x)
        depth = cur_y

        compose!(figs, fig)
        compose!(figs, label_fig)
    end 

    lines_ctx = compose!(context(order=-1), stroke("#333"))

    push!(queue, root)
    while !isempty(queue)
        cur = shift!(queue)
        pos = positions[cur]


        children = find(T[cur,:])

        for c in children
            push!(queue, c)

            c_pos = positions[c]
            compose!(lines_ctx, line([ pos, c_pos]))
        end


    end


    if Z != nothing
        zfigs = compose!(context(0,0), stroke("black"))

        for cur in keys(Z)
            label = labels == nothing ? "$cur" : labels[cur]

            fig = context(1, cur)
            compose!(fig, circle(0.5, 0.5, figsize/2), stroke("#333"), fill("#333"))
            label_fig = context(1, cur, order=1)
            compose!(label_fig, text(0.5, 0.5, label, hcenter, vcenter, Rotation()), stroke("white"), fill("white"), fontsize(8))

            compose!(zfigs, fig)
            compose!(zfigs, label_fig)

            text_fig = context(3, cur)
            compose!(text_fig, text(0.5, 0.5, Z[cur], hleft, vcenter, Rotation()), stroke("black"), fill("black"), fontsize(4))

            compose!(zfigs, text_fig)
        end
        tree = compose!(context(0,0,1,1,units=UnitBox(0, 0, max_x+1, depth+1)),
                   (context(order=-2), rectangle(), fill("white")),
                   lines_ctx, figs)

        panel = compose!(context(1,0,2,1,units=UnitBox(0,0, 30, length(Z)+2)),(context(order=-2), rectangle(), fill("white")), zfigs)

        full_fig = compose!(context(units=UnitBox(0,0,2,1)), tree, panel)

        return full_fig
    else
        tree = compose!(context(units=UnitBox(0, 0, max_x+1, depth+1)),
                   (context(order=-2), rectangle(), fill("white")),
                   lines_ctx, figs)
        return tree 
    end

end

function get_mutation_plot_locations(nU, x, y_min, y_max)
    x_inds = Float64[]
    y_inds = Float64[]

    dy = (y_max-y_min)/(nU+1)
    y_min += dy
    y_max -= dy 
    for y = linspace(y_min, y_max, nU)
        push!(x_inds, x)
        push!(y_inds, y)
    end
    (x_inds, y_inds)
end

function plot_Z_Y_pY(Z, Y, effects; plot=true, sortK=true)

    (N,K) = size(Z)

    if sortK
        mk = sum(Z,1)
        perm_K = sortperm(mk[:])
        Z = Z[:,perm_K]
    end
 
    Z, perm = uof(Z[:, end:-1:1])
    N = length(perm)

    perm = reverse(perm)
    Z = Z[reverse([1:N]),:]

    Y = deepcopy(Y)

    for s = 1:length(Y)
        YY = Y[s]
        YY[find(YY .< 0)] = 0
    end

    Y = mean(Y)
    Y = Y[perm,perm]
    pY = exp(broadcast(log_predictive, effects))
    pY = pY[perm,perm]

    colormap("grays")
    plot_Y = imagesc(Y[end:-1:1,:])
    plot_pY = imagesc(pY[end:-1:1,:])

    if K == 0
        Z = zeros(N,1)
    end

    plot_Z = imagesc(Z[end:-1:1,:])
    setattr(plot_Z,"title", "Z")
    setattr(plot_pY,"title", "p(Y)")
    setattr(plot_Y,"title", "Y")


    if plot
        t = Table(1,3)
        t[1,1] = plot_Z
        t[1,2] = plot_pY
        t[1,3] = plot_Y
        Winston.tk(t)
        t
    else
        plot_Z, plot_Y, plot_pY
    end
end

function plot_Z_Y_W(Z, Y, W; plot=true)

    (N,K) = size(Z)

    Z, perm = uof(Z)
    N = length(perm)

    perm = reverse(perm)
    Z = Z[reverse([1:N]),:]

    Y = deepcopy(Y)

    for s = 1:length(Y)
        YY = Y[s]
        YY[find(YY .< 0)] = 0
    end

    Y = mean(Y)
    Y = Y[perm,perm]
    Y = Y[end:-1:1,:]

    colormap("grays")
    plot_Y = imagesc(Y)
    setattr(plot_Y, "title", "Network Y")
    setattr(plot_Y.x1, "draw_ticklabels", false)
    setattr(plot_Y.y1, "draw_ticklabels", false)

    plot_W = imagesc(W[end:-1:1,:])
    setattr(plot_W, "title", "Weights W")
    setattr(plot_W, "xlabel", "Feature")
    setattr(plot_W, "ylabel", "Feature")
    setattr(plot_W.x1, "draw_ticklabels", false)
    setattr(plot_W.y1, "draw_ticklabels", false)

    if K == 0
        Z = zeros(N,1)
    end

    plot_Z = imagesc(Z[end:-1:1,:])
    setattr(plot_Z, "ylabel", "Actor")
    setattr(plot_Z, "title", "Features Z")
    setattr(plot_Z.x1, "draw_ticklabels", false)
    setattr(plot_Z.y1, "draw_ticklabels", false)

    plot_Zt = imagesc(Z'[end:-1:1,:])
    setattr(plot_Zt, "xlabel", "Actor")
    setattr(plot_Zt, "title", "Features Z^T")
    setattr(plot_Zt.x1, "draw_ticklabels", false)
    setattr(plot_Zt.y1, "draw_ticklabels", false)

    if plot
        t = Table(2,2)
        t[1,1] = plot_Z
        t[1,2] = plot_Y
        t[2,1] = plot_W
        t[2,2] = plot_Zt
        Winston.tk(t)
        t
    else
        plot_Z, plot_Y, plot_W, plot_Zt
    end
end

function createZfromTree(tree, features; respect_W=false)
    numActors = size(tree,1) + 1;
    Z = zeros(numActors, iround(sum(features)));

    # map from node to index in Z
    feature_indices = cumsum(features) - features + 1
    feature_indices[find(features == 0)] = 0

    Z = addFeaturesBelow(tree, features, Z, feature_indices, Int64[], 2 * numActors - 1);
    Z
end

function addFeaturesBelow(tree, features, Z, feature_indices, featuresAboveCurrentNode, currentNode)
    numActors = size(Z,1);

    for i = 1:features[currentNode]
        cur_feature = feature_indices[currentNode]+i-1
        featuresAboveCurrentNode = [featuresAboveCurrentNode, cur_feature];
    end
    
    if currentNode <= numActors 
        for i = 1:length(featuresAboveCurrentNode)
            Z[currentNode, featuresAboveCurrentNode[i]] = 1;
        end
        return Z 
    end

    Z = addFeaturesBelow(tree, features, Z, feature_indices, featuresAboveCurrentNode, iround(tree[currentNode-numActors, 1])); 
    Z = addFeaturesBelow(tree, features, Z, feature_indices, featuresAboveCurrentNode, iround(tree[currentNode-numActors, 2])); 
    
    Z 
end

function make_monk_labels()

    labels = cell(18);
    labels[1] = "w";
    labels[2] = "o";
    labels[3] = "o";
    labels[4] = "o";
    labels[5] = "o";
    labels[6] = "o";
    labels[7] = "w";

    labels[8] = "y";
    labels[9] = "y";
    labels[10] = "y";
    labels[11] = "y";
    labels[12] = "y";
    labels[13] = "y";
    labels[14] = "y";
    labels[15] = "w";

    labels[16] = "+";
    labels[17] = "+";
    labels[18] = "+";
    labels
end

function box_plot_along_dim(data::Array{Float64}, dim::Int64; offset=0.0, p = FramedPlot())

    nd = ndims(data)

    N = size(data,dim)

    labels = [string(i) for i = 1:N]
    dict = Dict{ASCIIString, Vector{Float64}}()

    for i = 1:N

        indices = [ [1:size(data,j)] for j = 1:nd]
        indices[dim] = [i]

        dict[labels[i]] = data[indices...][:]
    end

    p = box_plot(dict, order=labels, labels=labels, offset=offset, p = p)
end

function box_plot(data::Dict{ASCIIString, Vector{Float64}};
                  labels::Vector{ASCIIString} = [k for k in keys(data)],
                  order::Vector{ASCIIString} = [k for k in keys(data)],
                  offset::Float64=0.0,
                  p = FramedPlot())
    vals::Vector{Vector{Float64}} = [data[k] for k in order]
    p = box_plot(vals, offset, p=p)
    setattr(p.x1,"ticklabels", labels)
    p
end

function box_plot(data::Vector{Vector{Float64}}, offset::Float64; p = FramedPlot())

    width = 0.1
    cl = 0x3083FF
    for i = 1:length(data)
        Q3 = quantile(data[i],0.75)
        Q2 = quantile(data[i],0.5)
        Q1 = quantile(data[i],0.25)

        IQR = Q3-Q1

        lower_whisker = Q1-1.5*IQR
        upper_whisker = Q3+1.5*IQR

        u = find(data[i] .< upper_whisker)
        upper_whisker = length(u) > 0 ? maximum(data[i][u]) : maximum(data[i])

        l = find(data[i] .> lower_whisker)
        lower_whisker = length(l) > 0 ? minimum(data[i][l]) : minimum(data[i])

        X = i+offset
        Xl = X-width
        Xr = X+width

        whisker_mean = (upper_whisker+lower_whisker)/2.0
        whisker_err = upper_whisker-whisker_mean
        add(p, SymmetricErrorBarsY(X, whisker_mean, whisker_err))
        add(p, FillBetween( [Xl, Xr], [Q1, Q1], [Xl, Xr], [Q3, Q3], "color", cl ))
        add(p, Curve([Xl, Xr] , [Q2, Q2]))

        outliers = find( data[i] .< lower_whisker)
        outliers = [outliers, find(data[i] .> upper_whisker)]
        if length(outliers) > 0
            add(p, Points(X*ones(length(outliers)), data[i][outliers], "symboltype", "asterisk"))
        end
    end
    setattr(p.x1, "ticks", [1.0:length(data)])
    p
end

