using Winston
import Winston.Interval
function imagegs{T<:Real}(xrange::Interval, yrange::Interval, data::AbstractArray{T,2}, clims::Interval; plot=true)
    p = FramedPlot()
    setattr(p, "xrange", xrange)
    setattr(p, "yrange", reverse(yrange))
    img = Winston.data2rgb(data, clims, Winston.GrayColormap())
    add(p, Image(xrange, reverse(yrange), img))

    if plot
        Winston.tk(p)
    else
        return p
    end
end

imagegs(xrange, yrange, data; plot=true) = imagegs(xrange, yrange, data, (min(data),max(data)), plot=plot)
imagegs(data; plot=true) = ((h, w) = size(data); imagegs((0,w), (0,h), max(data)-data, plot=plot))

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

function dendrogram(Z, U; plot=true)

    Nm1, _ = size(Z)
    N = Nm1+1
    # Sort in depth-first-search order so the dendrogram doesn't cross itself 
    times = zeros(2N-1)
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
        
        push!(visited_stack, l)
        push!(visited_stack, r) 
    end

    locations[leaf_order] = [1:N]/N

    # Traverse in breadth first order so get leaves to root ordering of internal nodes
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
                n_x, n_y = get_mutation_plot_locations(U[n], locations[n], n_t, 1.0)
                append!(mutations_x, n_x)
                append!(mutations_y, n_y)
            end
        end
    end

    add(p, Points(mutations_x, mutations_y, "type", "diamond"))

    if plot
        Winston.tk(p)
    end
    p
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

function plot_Z_Y_pY(Z, Y, effects; plot=true)

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
    pY = broadcast(log_predictive, effects)
    pY = pY[perm,perm]

    plot_Y = imagegs(Y, plot=false)
    plot_Z = imagegs(Z, plot=false)
    plot_pY = imagegs(pY, plot=false)


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
