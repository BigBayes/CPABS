abstract Node

type Nil <: Node
end

type TreeNode{T <: Real} <: Node
    location::Array{T,1}
    parent::Node
    children::Array{Node,1}
    index::Int64
    num_leaves::Int64 #more accurately number of leaves
    num_ancestors::Int64 #not self inclusive
end

function TreeNode{T}(loc::Array{T,1}, ind::Int64)
    TreeNode{T}(loc,Nil(),Array(Node,2),ind,1,0)
end

show(tree::TreeNode) = print(tree.index)

type Tree{T}
    nodes::Array{TreeNode{T},1}
    Tree() = new()
end


function Tree(U::Array{Int64,1})
    _2nm1 = size(U)
    N = (_2nm1+1)/2
    tree = Tree{Int64}()
    tree.nodes = Array(TreeNode{Int64},2N-1)

    for i = 1:N
        tree.nodes[i] = TreeNode(U[i],i)
        tree.nodes[i].children[1] = Nil();
        tree.nodes[i].children[2] = Nil();
    end
    coalescing_nodes = linspace(1,N,N)
    coalescing_nodes_remaining = N

    for i = N+1:2N-1
        tree.nodes[i] = TreeNode(U[i], i)

        #coalesce pairs uniformly at random
        l_ind = randi(coalescing_nodes_remaining)
        l = coalescing_nodes[l_ind]
        r_ind = l_ind
        while r_ind == l_ind
            r_ind = randi(coalescing_nodes_remaining)
        end
        r = coalescing_nodes[r_ind]

        #update the set of nodes remaining
        m_ind = min(l_ind,r_ind)
        x_ind = max(l_ind,r_ind)
        del(coalescing_nodes, x_ind)
        del(coalescing_nodes, m_ind)
        push(coalescing_nodes, i)
        coalescing_nodes_remaining -= 1
     
        tree.nodes[l].parent = tree.nodes[i]
        tree.nodes[r].parent = tree.nodes[i]

        tree.nodes[i].children = [tree.nodes[l], tree.nodes[r]]

        tree.nodes[i].num_leaves = tree.nodes[l].num_leaves +
                                        tree.nodes[r].num_leaves 
    end

    for i = 2N-2:-1:1
        tree.nodes[i].num_ancestors = tree.nodes[i].parent.num_ancestors + 1
    end 

    tree 
end


function PruneIndexFromTree!{T}(tree::Tree{T}, index::Int)
    self = tree.nodes[index]
    parent = self.parent
    grandparent = parent.parent

    assert(grandparent != Nil())

    #is the parent the left or right child of grandparent
    parent_direction = find( grandparent.children .== parent )[1]
    parent_sibling = 
        grandparent.children[find( grandparent.children .!= parent)[1]]

    self_direction = find( parent.children .== self)[1]
    sibling_direction = find( parent.children .!= self)[1]
    sibling = parent.children[sibling_direction]

    grandparent.children[parent_direction] = sibling
    sibling.parent = grandparent

    parent.parent = Nil()
    parent.children[sibling_direction] = Nil()

    UpdateDescendantCounts(tree, grandparent)

    UpdateSubtreeAncestorCounts(tree, sibling)
    UpdateSubtreeAncestorCounts(tree, parent)

    nothing 
end

#insert the subtree rooted at insert_index above the node at location_index
function InsertIndexIntoTree!{T}(tree::Tree{T}, 
                                 insert_index::Int,
                                 location_index::Int)
    self = tree.nodes[insert_index]
    parent = self.parent
    new_sibling = tree.nodes[location_index]
    grandparent = new_sibling.parent

    self_direction = find(parent.children .== self)[1]
    sibling_direction = find(parent.children .== Nil())[1]

    new_sibling_direction = find(grandparent.children .== new_sibling)[1]

    parent.children[sibling_direction] = new_sibling 
    grandparent.children[new_sibling_direction] = parent

    new_sibling.parent = parent
    parent.parent = grandparent

    UpdateDescendantCounts(tree, parent)
    UpdateSubtreeAncestorCounts(tree, parent)
    nothing
end

function UpdateDescendantCounts!{T}(tree::Tree{T},
                                    start_node::TreeNode{T})
    cur = start_node
    while cur != Nil()
        cur.num_leaves = cur.children[1].num_leaves 
                            + cur.children[2].num_leaves 
        cur = cur.parent
    end
end

function UpdateSubtreeAncestorCounts!{T}(tree::Tree{T},
                                         subtree_root::TreeNode{T})
    subtree_queue = TreeNode{T}[]
    enqueue(subtree_queue, subtree_root)
    while length(subtree_queue) > 0
        cur = pop(subtree_queue)
        if cur != Nil() 
            enqueue(subtree_queue, cur.children[1])
            enqueue(subtree_queue, cur.children[2])
            cur.num_ancestors = cur.parent.num_ancestors + 1
        end
    end
end

# Returns a vector of nodes indicating the path from leaf to root
function GetPath{T}(tree::Tree{T},
                    node_index::Int)
    path = Array(TreeNode{T},0)
    current_node = tree.nodes[node_index]
    push(path, current_node)
    while (current_node != Nil())
        push(path, current_node)
        current_node = current_node.parent
    end
    path
end

function GetPathNumLeaves{T}(tree::Tree{T},
                             node_index::Int)
    path_size = 0
    current_node = tree.nodes[node_index]
    while (current_node != Nil())
        path_size += 1 
        current_node = current_node.parent
    end

    num_leaves = zeros(Int, path_size)

    current_node = tree.nodes[node_index]
    path_index = 1
    while (current_node != Nil())
        num_leaves[path_index] = current_node.num_leaves
        current_node = current_node.parent
        path_index += 1
    end

    num_leaves
end

function ProposedNumLeaves{T}(tree::Tree{T},
                              pruned_index::Int,
                              insert_index::Int)
    current_node = tree.nodes[insert_index]
    num_leaves_added = tree.nodes[pruned_index].num_leaves

    path_num_leaves = GetPathNumLeaves(tree,insert_index)
    
    path_num_leaves + num_leaves_added
end

function GetLeaves{T}(tree::Tree{T},
                      subtree_root::Int)

    if index <= (length(tree.nodes)+1)/2
        return [index]
    else
        a = Array(Int64,0)
        left = tree.nodes[index].children[1]
        right = tree.nodes[index].children[2]
        if left != Nil()
            append!(a,GetLeaves(tree,left.index))
        end
        if right != Nil()
            append!(a,GetLeaves(tree,right.index))
        end
        return a
    end

end

#don't think this works right now
function GetRandomLeaf{T}(tree::Tree{T},
                                index::Int)
    leaves = GetLeaves(tree,index)
    d_ind = randi(length(leaves))
    leaves[d_ind]
end

function FindRoot{T}(tree::Tree{T},
                     start_index::Int64)
    cur = tree.nodes[start_index]
    while cur.parent != Nil()
        cur = cur.parent
    end
    cur
end

function GetSubtreeIndicies{T}(tree::Tree{T},
                               subtree_root::Int64)
    subtree_indices = IntSet()
    queue = Int64[]
    enqueue(queue, subtree_root)
    while length(queue) > 0
        cur = tree.nodes[pop(queue)]
        if cur != Nil()
            enqueue(queue, cur.children[1].index)
            enqueue(queue, cur.children[2].index)
            add(subtree_indices, cur.index)
        end  
    end

    subtree_indices
end

# stack will contain a leaves to root ordering of nodes in the pruned tree
function GetLeafToRootOrdering{T}(tree::Tree{T},
                                  root_index::Int64)
    queue = Int64[]
    stack = Int64[]
    enqueue(queue, root_index)

    while length(queue) > 0
        cur = tree.nodes[pop(queue)]
        if cur != Nil()
            enqueue(queue, cur.children[1].index)
            enqueue(queue, cur.children[2].index)
            push(stack, cur.index)
        end  
    end

    reverse(stack)
end

function ConstructZ{T}(tree::Tree{T})
    N = (length(tree.nodes) + 1) / 2
    U = zeros(Int64, 2N - 1)
    for i = 1:2N-1
        U[i] = tree.nodes[i].location
    end

    Z = zeros(Int64, N, sum(U))

    cur_j = 0
    for i = 1:2N-1
        leaves = GetLeaves(tree, i)
        for j = 1:U[i]
            cur_j += 1
            Z[leaves, cur_j] = 1
        end
    end
    Z
end

