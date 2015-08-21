import Base.copy
import Base.show
import Base.==
import Base.serialize
abstract Node

type Nil <: Node
end

type TreeNode{T} <: Node
    state::T
    parent::Node
    children::Array{Node,1}
    index::Int64
    num_leaves::Int64 #number of leaves in subtree
    num_ancestors::Int64 #self inclusive!
    rhot::Float64 #proportion of immediate parent in *infinite* tree
                 #shared with sibling
    rho::Float64 #splitting point in Aldous beta-splitting representation, 
                 #represented as the prob. mass of the left branch conditioned on rhot
end

function TreeNode{T}(state::T, ind::Int64)
    TreeNode{T}(state,Nil(),Array(Node,2),ind,1,0,1.0,0.5)
end

show(tree_node::TreeNode) = print(tree_node.index)
copy(tree_node::TreeNode) = TreeNode(deepcopy(tree_node.state), tree_node.index)
function ==(n1::TreeNode, n2::TreeNode)
    b1 = n1.state == n2.state && n1.index == n2.index
    if n1.parent == Nil() || n2.parent == Nil()
        b2 = n1.parent == n2.parent
    else
        b2 = n1.parent.index == n2.parent.index
    end

    if n1.children[1] == Nil() || n2.children[1] == Nil()
        b3 = n1.children[1] == n2.children[1] && n2.children[2] == n2.children[2]
    else
        b3 = n1.children[1].index == n2.children[1].index &&
             n2.children[2].index == n2.children[2].index
    end

    b1 && b2 && b3
end

type Tree{T}
    nodes::Array{TreeNode{T},1}
    Tree() = new()

end

function Tree(t::Array{TreeNode{Int64},1})
    tree = Tree{Int64}()
    tree.nodes = t
    tree
end

function ==(tree1::Tree, tree2::Tree)
    b = length(tree1.nodes) == length(tree2.nodes)
    if !b
        return b
    end

    for i = 1:length(tree1.nodes)
        b = b && tree1.nodes[i] == tree2.nodes[i]
    end
    b
end

function copy(tree::Tree{Int64})
    n = length(tree.nodes)
    new_tree = Tree{Int64}()
    new_tree.nodes = Array(TreeNode{Int64}, n)  
    for i = 1:n
        new_tree.nodes[i] = copy(tree.nodes[i])
        new_tree.nodes[i].num_ancestors = tree.nodes[i].num_ancestors
        new_tree.nodes[i].num_leaves = tree.nodes[i].num_leaves
        new_tree.nodes[i].rho = tree.nodes[i].rho
        new_tree.nodes[i].rhot = tree.nodes[i].rhot
    end

    for i = 1:n
        if tree.nodes[i].parent == Nil()
            p = Nil()
        else
            p = new_tree.nodes[tree.nodes[i].parent.index]
        end

        if tree.nodes[i].children[1] == Nil()
            l = Nil()
            r = Nil()
        else
            l = new_tree.nodes[tree.nodes[i].children[1].index]
            r = new_tree.nodes[tree.nodes[i].children[2].index]
        end
        new_tree.nodes[i].parent = p
        new_tree.nodes[i].children = [l, r]
    end
    new_tree
end

function Tree(U::Array{Int64,1})
    (_2nm1,) = size(U)
    N::Int64 = (_2nm1+1)/2
    tree = Tree{Int64}()
    tree.nodes = Array(TreeNode{Int64},2N-1)

    for i = 1:N
        tree.nodes[i] = TreeNode(U[i],i)
        tree.nodes[i].children[1] = Nil();
        tree.nodes[i].children[2] = Nil();
    end
    coalescing_nodes = [1:N]
    coalescing_nodes_remaining = N

    for i = N+1:2N-1
        tree.nodes[i] = TreeNode(U[i], i)

        #coalesce pairs uniformly at random
        l_ind = rand(1:coalescing_nodes_remaining)
        l = coalescing_nodes[l_ind]
        r_ind = l_ind
        while r_ind == l_ind
            r_ind = rand(1:coalescing_nodes_remaining)
        end
        r = coalescing_nodes[r_ind]

        #update the set of nodes remaining
        m_ind = min(l_ind,r_ind)
        x_ind = max(l_ind,r_ind)
        splice!(coalescing_nodes, x_ind)
        splice!(coalescing_nodes, m_ind)
        push!(coalescing_nodes, i)
        coalescing_nodes_remaining -= 1
   
#        println("i,l,r: $i $l $r")
#        println("tree.nodes[i]: $(tree.nodes[i])") 
#        println("tree.nodes[l]: $(tree.nodes[l])") 
#        println("tree.nodes[r]: $(tree.nodes[r])") 
#
 
        tree.nodes[l].parent = tree.nodes[i]
        tree.nodes[r].parent = tree.nodes[i]

        tree.nodes[i].children = [tree.nodes[l], tree.nodes[r]]

        tree.nodes[i].num_leaves = tree.nodes[l].num_leaves +
                                   tree.nodes[r].num_leaves 
    end


    tree.nodes[2N-1].num_ancestors = 1
    for i = 2N-2:-1:1
        tree.nodes[i].num_ancestors = tree.nodes[i].parent.num_ancestors + 1
    end 

    tree 
end

function Tree(U::Array{Vector{Float64},1})
    (_2nm1,) = size(U)
    N::Int64 = (_2nm1+1)/2
    tree = Tree{Vector{Float64}}()
    tree.nodes = Array(TreeNode{Vector{Float64}},2N-1)

    for i = 1:N
        tree.nodes[i] = TreeNode(U[i],i)
        tree.nodes[i].children[1] = Nil();
        tree.nodes[i].children[2] = Nil();
    end
    coalescing_nodes = [2:N]
    coalescing_nodes_remaining = N-1

    for i = N+1:2N-2
        tree.nodes[i] = TreeNode(U[i], i)

        #coalesce pairs uniformly at random
        l_ind = rand(1:coalescing_nodes_remaining)
        l = coalescing_nodes[l_ind]
        r_ind = l_ind
        while r_ind == l_ind
            r_ind = rand(1:coalescing_nodes_remaining)
        end
        r = coalescing_nodes[r_ind]

        #update the set of nodes remaining
        m_ind = min(l_ind,r_ind)
        x_ind = max(l_ind,r_ind)
        splice!(coalescing_nodes, x_ind)
        splice!(coalescing_nodes, m_ind)
        push!(coalescing_nodes, i)
        coalescing_nodes_remaining -= 1
   
#        println("i,l,r: $i $l $r")
#        println("tree.nodes[i]: $(tree.nodes[i])") 
#        println("tree.nodes[l]: $(tree.nodes[l])") 
#        println("tree.nodes[r]: $(tree.nodes[r])") 
#
 
        tree.nodes[l].parent = tree.nodes[i]
        tree.nodes[r].parent = tree.nodes[i]

        tree.nodes[i].children = [tree.nodes[l], tree.nodes[r]]

        tree.nodes[i].num_leaves = tree.nodes[l].num_leaves +
                                   tree.nodes[r].num_leaves 
    end

    # We require that the left subtree from the root is a leaf 
    # for the purposes of the subclonal reconstruction model, we don't expect two independent cancer populations, etc

    tree.nodes[2N-1] = TreeNode(U[2N-1], 2N-1)
    l = 1
    r = 2N-2
    tree.nodes[l].parent = tree.nodes[2N-1]
    tree.nodes[r].parent = tree.nodes[2N-1]

    # likelihood functions assume right child is first, and here it matters 
    tree.nodes[2N-1].children = [tree.nodes[r], tree.nodes[l]]

    tree.nodes[2N-1].num_leaves = tree.nodes[l].num_leaves +
                                  tree.nodes[r].num_leaves 
    


    tree.nodes[2N-1].num_ancestors = 1
    for i = 2N-2:-1:1
        tree.nodes[i].num_ancestors = tree.nodes[i].parent.num_ancestors + 1
    end 

    tree 
end

function Tree(U::Array{Vector{Float64},1},
              coalescing_order::Array{Int64,2})
    (_2nm1,) = size(U)
    N::Int64 = (_2nm1+1)/2
    tree = Tree{Vector{Float64}}()
    tree.nodes = Array(TreeNode{Vector{Float64}},2N-1)

    for i = 1:N
        tree.nodes[i] = TreeNode(U[i],i)
        tree.nodes[i].children[1] = Nil();
        tree.nodes[i].children[2] = Nil();
    end

    for i = N+1:2N-2
        tree.nodes[i] = TreeNode(U[i], i)

        l = coalescing_order[i-N,1]
        r = coalescing_order[i-N,2]

 
        tree.nodes[l].parent = tree.nodes[i]
        tree.nodes[r].parent = tree.nodes[i]

        tree.nodes[i].children = [tree.nodes[l], tree.nodes[r]]

        tree.nodes[i].num_leaves = tree.nodes[l].num_leaves +
                                   tree.nodes[r].num_leaves 
    end

    # We require that the left subtree from the root is a leaf 
    # for the purposes of the subclonal reconstruction model, we don't expect two independent cancer populations, etc

    tree.nodes[2N-1] = TreeNode(U[2N-1], 2N-1)
    l = 1
    r = 2N-2
    tree.nodes[l].parent = tree.nodes[2N-1]
    tree.nodes[r].parent = tree.nodes[2N-1]

    # likelihood functions assume right child is first, and here it matters 
    tree.nodes[2N-1].children = [tree.nodes[r], tree.nodes[l]]

    tree.nodes[2N-1].num_leaves = tree.nodes[l].num_leaves +
                                  tree.nodes[r].num_leaves 
    


    tree.nodes[2N-1].num_ancestors = 1
    for i = 2N-2:-1:1
        tree.nodes[i].num_ancestors = tree.nodes[i].parent.num_ancestors + 1
    end 

    tree 
end


function copy(tree::Tree{Vector{Float64}})
    n = length(tree.nodes)
    new_tree = Tree{Vector{Float64}}()
    new_tree.nodes = Array(TreeNode{Vector{Float64}}, n)  
    for i = 1:n
        new_tree.nodes[i] = copy(tree.nodes[i])
        new_tree.nodes[i].num_ancestors = tree.nodes[i].num_ancestors
        new_tree.nodes[i].num_leaves = tree.nodes[i].num_leaves
        new_tree.nodes[i].rho = tree.nodes[i].rho
        new_tree.nodes[i].rhot = tree.nodes[i].rhot
    end

    for i = 1:n
        if tree.nodes[i].parent == Nil()
            p = Nil()
        else
            p = new_tree.nodes[tree.nodes[i].parent.index]
        end

        if tree.nodes[i].children[1] == Nil()
            l = Nil()
            r = Nil()
        else
            l = new_tree.nodes[tree.nodes[i].children[1].index]
            r = new_tree.nodes[tree.nodes[i].children[2].index]
        end
        new_tree.nodes[i].parent = p
        new_tree.nodes[i].children = [l, r]
    end
    new_tree
end

function SwapIndices!{T}(tree::Tree{T}, index1::Int64, index2::Int64)
    node1 = tree.nodes[index1]
    node2 = tree.nodes[index2]

    node1.index = index2
    node2.index = index1

    tree.nodes[index1] = node2
    tree.nodes[index2] = node1

end

# Construct a new tree object from a specifed node of an existing tree 
# (eg after pruning or to obtain a subtree)
# Can also be used to reindex the nodes (eg after grafting a new subtree in),
# as long as all node indices are unique (the leaves don't need to have the first N indices)
function MakeReindexedTree{T}(old_tree::Tree{T}, root_index::Int64, old_N::Int64)
    tree = Tree{T}()

    old_root = FindRoot(old_tree, root_index)
    old_leaves = GetLeaves(old_tree, old_root.index)


    root_node = old_tree.nodes[root_index]
    leaves = GetLeaves(old_tree, root_index)
    N = length(leaves)
    tree.nodes = Array(TreeNode{T}, 2N-1)

    queue = copy(leaves)

    index_map = zeros(Int64, 2max(old_N,N)-1)

    # build index_map and instantiate nodes of the new tree
    println("N: $N")
    indices = [1:2N-1]
    while length(queue) > 0
        println("queue: $queue")
        node_index = shift!(queue)

        println("node_index: $node_index")
        println("indices: $indices")
        new_index = shift!(indices)

 
        node = old_tree.nodes[node_index]
        new_node = copy(node)
        new_node.rhot = node.rhot 
        new_node.rho = node.rho 
        new_node.children[1] = Nil()
        new_node.children[2] = Nil()

        tree.nodes[new_index] = new_node

        if node.parent != Nil() && node_index != root_index
            parent_index = node.parent.index
            if index_map[parent_index] == 0 && !(parent_index in queue)
                push!(queue, parent_index)
            end
        end

    end


    stack = Int64[]
    push!(stack, root_index)

    # connect parents/children in the new tree
    while length(stack) > 0
        node_index = pop!(stack)
        node = old_tree.nodes[node_index]

        new_index = index_map[node_index]
        new_node = tree.nodes[new_index]
        new_node.index = new_index


        l = node.children[2]
        r = node.children[1]
        if l != Nil()
            l_index = index_map[l.index]            
            new_node.children[2] = tree.nodes[l_index]
            tree.nodes[l_index].parent = new_node

            push!(stack, l.index)
        end
        if r != Nil()
            r_index = index_map[r.index]            
            new_node.children[1] = tree.nodes[r_index]
            tree.nodes[r_index].parent = new_node

            push!(stack, r.index)
        end

    end 

    root = FindRoot(tree, 1)  

    @assert root.index == index_map[root_index]
    println("root.index: $(root.index)")
    println("N: $N")
    UpdateSubtreeDescendantCounts!(tree, tree.nodes[index_map[root_index]])
    UpdateSubtreeAncestorCounts!(tree, tree.nodes[index_map[root_index]])
    return (tree, index_map[1:2*old_N-1])
end

function PruneIndexFromTree!{T}(tree::Tree{T}, index::Int)
    self = tree.nodes[index]
    parent = self.parent
    grandparent = parent.parent

    if grandparent != Nil()
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

        UpdateDescendantCounts!(tree, grandparent)

        UpdateSubtreeAncestorCounts!(tree, sibling)
        UpdateSubtreeAncestorCounts!(tree, parent)
    else
        self_direction = find( parent.children .== self)[1]
        sibling_direction = find( parent.children .!= self)[1]
        sibling = parent.children[sibling_direction]

        sibling.parent = Nil()
        if sibling.children[1] != Nil() 
            UpdateDescendantCounts!(tree, sibling)
        end
        parent.children[sibling_direction] = Nil()

        UpdateSubtreeAncestorCounts!(tree, sibling)
        UpdateSubtreeAncestorCounts!(tree, parent)
    end
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

    if all(parent.children .!= Nil())
        println("child1: $(parent.children[1].index)")
        println("child2: $(parent.children[2].index)")
    end
    sibling_direction = find(parent.children .== Nil())[1]

    parent.children[sibling_direction] = new_sibling 
    if grandparent != Nil()
        new_sibling_direction = find(grandparent.children .== new_sibling)[1]
        grandparent.children[new_sibling_direction] = parent
    end

    new_sibling.parent = parent
    parent.parent = grandparent

    UpdateDescendantCounts!(tree, parent)
    UpdateSubtreeAncestorCounts!(tree, parent)
    nothing
end

# Remove index from tree, promoting the right child
function RemoveAndPromoteRight!{T}(tree::Tree{T},
                                   removal_index::Int)

    self = tree.nodes[removal_index]
    parent = self.parent


    right_child = self.children[1]
    left_child = self.children[2]

    if parent != Nil()
        self_direction = find(parent.children .== self)[1]
        parent.children[self_direction] = right_child
        right_child.parent = parent
    else
        right_child.parent = Nil()
    end

    right_left_child = right_child.children[2]

    rl_cur = right_left_child
    while rl_cur != Nil() && rl_cur.children[2] != Nil()
        rl_cur = rl_cur.children[2]
    end 

    if rl_cur != Nil()
        rl_cur = rl_cur.parent
        rl_cur.children[2].parent = Nil() #remove this leaf node
        rl_cur.children[2] = left_child
        left_child.parent = rl_cur    
    else
        if parent == Nil()
            left_child.parent = Nil()
        else
            parent.children[self_direction] = left_child
            left_child.parent = parent
        end
    end

    UpdateDescendantCounts!(tree, left_child.parent)
    UpdateSubtreeAncestorCounts!(tree, left_child.parent)
    nothing
end

# As there are multiples trees after which using RemoveAndPromoteRight! will
# result in the same tree, there are multiple ways in which we may build a 
# tree by inserting a node at a particular location.  
# That is to remove node "n," and promote "r," we need to move node "rl" above "l"
#            
#
#          ...                    ...          
#             \                      \
#              n                      r
#            /   \       =>         /   \
#          l       r              rl     o
#         / \     / \            / \    
#        o   o   rl  o          l   o   
#               / \            / \  
#              o   o          o   o
#
# where one of the leaves has been deleted
# However, we may have arrived at the same tree by removing "n" from these trees
#                                     ...
#                                        \
#            ...                          n          
#               \                        / \
#                n                      o   r
#              /   \                      /   \
#            rl     r                   rl     o
#           / \    / \                 / \    
#          l   o  o   o               l   o   
#         / \                        / \  
#        o   o                      o   o
#        
# So in order to map back to all possible trees that can result in the current tree by removing
# the newly inserted node, we have to MOVE some subset of the internal nodes on the path to the
# leftmost leaf of the right child of the newly inserted node TO the left subtree of the newly 
# inserted node. N_move determines the number that are moved.   
#
# Assumes that tree.nodes[insert_index] has a left child that is not Nil() 
function InsertAndDemoteLeft!{T}(tree::Tree{T},
                                 insert_index::Int,
                                 location_index::Int,
                                 N_move::Int)

    self = tree.nodes[insert_index]
    left_child = self.children[2]
    
    @assert left_child != Nil()
 
    right_child = tree.nodes[location_index]
    right_left_child = right_child.children[2]
    rl_prune = right_left_child

    if N_move > 0
        N_count = N_move
        while N_count > 0
            N_count -= 1
            rl_prune = rl_prune.children[2]        
        end

        @assert rl_prune != Nil()

        # remove the relevent subtree
        rl_prune_child = rl_prune.children[2]
        rl_prune_child.parent = right_child
        right_child.children[2] = rl_prune_child

        # attach it to the left of self
        left_child.parent = rl_prune
        rl_prune.children[2] = left_child
        right_left_child.parent = self
        self.children[2] = right_left_child
    end

    self.children[1] = right_child
    self.parent = right_child.parent
    right_child.parent = self 
    nothing
end

function UpdateDescendantCounts!{T}(tree::Tree{T},
                                    start_node::TreeNode{T})
    cur = start_node
    while cur != Nil()
        cur.num_leaves = cur.children[1].num_leaves + 
                         cur.children[2].num_leaves 
        cur = cur.parent
    end
end

function UpdateSubtreeDescendantCounts!{T}(tree::Tree{T},
                                           subtree_root::TreeNode{T})
    if subtree_root == Nil()
        return
    end

    indices = GetLeafToRootOrdering(tree, subtree_root.index)

    for index in indices
        cur = tree.nodes[index]
        cur.num_leaves = 0
        l = cur.children[1]
        r = cur.children[2]
        if l == Nil() && r == Nil()
            cur.num_leaves = 1
        end
        if l != Nil()
            cur.num_leaves += l.num_leaves
        end
        if r != Nil()
            cur.num_leaves += r.num_leaves
        end
    end
end

function UpdateSubtreeAncestorCounts!{T}(tree::Tree{T},
                                         subtree_root::TreeNode{T})
    subtree_queue = Int64[]
    if subtree_root == Nil()
        return
    end

    unshift!(subtree_queue, subtree_root.index)
    while length(subtree_queue) > 0
        cur = tree.nodes[pop!(subtree_queue)]
        if cur.children[1] != Nil() 
            unshift!(subtree_queue, cur.children[1].index)
        end
        if cur.children[2] != Nil()
            unshift!(subtree_queue, cur.children[2].index)
        end
        if cur.parent == Nil()
            cur.num_ancestors = 1
        else
            cur.num_ancestors = cur.parent.num_ancestors + 1
        end
    end
end

# Returns a vector of nodes indicating the path from leaf to root
function GetPath{T}(tree::Tree{T},
                    node_index::Int)
    path = Array(Int64,0)
    current_node = tree.nodes[node_index]
    while (current_node != Nil())
        push!(path, current_node.index)
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

    if tree.nodes[subtree_root].children[1] == Nil() &&
       tree.nodes[subtree_root].children[2] == Nil()
        return [subtree_root]
    else
        a = Array(Int64,0)
        left = tree.nodes[subtree_root].children[1]
        right = tree.nodes[subtree_root].children[2]
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
    d_ind = rand(1:length(leaves))
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
    assert( tree.nodes[subtree_root] != Nil())
    queue = Int64[]
    unshift!(queue, subtree_root)
    while length(queue) > 0
        cur = tree.nodes[pop!(queue)]
        if cur.children[1] != Nil()
            unshift!(queue, cur.children[1].index)
        end
        if cur.children[2] != Nil()
            unshift!(queue, cur.children[2].index)
        end
        push!(subtree_indices, cur.index)
    end

    subtree_indices
end

# Self inclusive
function GetAncestors{T}(tree::Tree{T},
                         node_index::Int64)
    ancestors = Array(TreeNode, 0)

    cur = tree.nodes[node_index]

    while cur != Nil()
        push!(ancestors, cur)
        cur = cur.parent
    end

    return ancestors
end

function GetDescendants{T}(tree::Tree{T},
                           node_index::Int64)
    descendants = Array(TreeNode, 0)
    stack = Array(TreeNode, 0)

    cur = tree.nodes[node_index]
    push!(stack, cur)

    while length(stack) > 0 
        c = pop!(stack)
        push!(descendants, c)

        r = c.children[1]
        l = c.children[2]

        if r != Nil()
            push!(stack, r)
        end

        if l != Nil()
            push!(stack, l)
        end
    end

    return descendants
end
# stack will contain a leaves to root ordering of nodes in the pruned tree
function GetLeafToRootOrdering{T}(tree::Tree{T},
                                  root_index::Int64)
    queue = Int64[]
    stack = Int64[]
    assert( tree.nodes[root_index] != Nil())
    unshift!(queue, root_index)

    while length(queue) > 0
        cur = tree.nodes[pop!(queue)]
        if cur.children[1] != Nil()
            unshift!(queue, cur.children[1].index)
        end
        if cur.children[2] != Nil()
            unshift!(queue, cur.children[2].index)
        end  
        push!(stack, cur.index)
    end

    reverse(stack)
end

function InitializeBetaSplits{T}(tree::Tree{T},
                                 draw_split::Function)

    N::Int = (length(tree.nodes) + 1) / 2

    for i = N+1:2N-1
        tree.nodes[i].rho = draw_split()
        tree.nodes[i].rhot = 1.0
    end
end

function UpdateBetaSplits{T}(tree::Tree{T},
                             update_split::Function)
    N::Int = (length(tree.nodes) + 1) / 2

    for i = N+1:2N-1
        cur = tree.nodes[i]
        left_child = cur.children[1]
        right_child = cur.children[2]

        cur.rho = update_split(left_child.num_leaves, right_child.num_leaves)
    end
end

function ConstructZ{T}(tree::Tree{T})
    N::Int = (length(tree.nodes) + 1) / 2
    U = zeros(Int64, 2N - 1)
    for i = 1:2N-1
        U[i] = tree.nodes[i].state
    end

    Z_I = Int64[]
    Z_J = Int64[]
    Z_S = Int64[]

    cur_j = 0
    for i = 1:2N-1
        leaves = GetLeaves(tree, i)
        for j = 1:U[i]
            cur_j += 1
            for l = leaves
                push!(Z_I, l)
                push!(Z_J, cur_j)
                push!(Z_S, 1)
            end
        end
    end
    full(sparse(Z_I, Z_J, Z_S, N, sum(U)))
end

# Is a an ancestor of b from the right subtree?
function IsRightAncestor{T}(tree::Tree{T},
                            index_a::Int64,
                            index_b::Int64)

    node_a = tree.nodes[index_a]
    node_b = tree.nodes[index_b]

    if node_a == node_b
        return true
    end

    last_node = node_b

    while node_b != Nil()
        if node_b == node_a
            if last_node == node_a.children[1]
                return true
            else
                return false
            end
        end
        last_node = node_b
        node_b = node_b.parent
    end

    return false
end

# Get adjacency matrix for the "right descendents are children"
# interpretation of the tree
function GetAdjacencyMatrix{T}(tree::Tree{T})
    N::Int = (length(tree.nodes) + 1) / 2
   
    M = zeros(N-1, N-1)
 
    for i = 1:N
        cur = tree.nodes[i]

        cur_direction = 0

        indices = Int64[]

        while cur_direction != 1 && cur.parent != Nil()
            par = cur.parent
            cur_direction = find(par.children .== cur)[1]

            if cur_direction == 1
                M[par.index-N, indices] = 1
            elseif cur_direction == 2
                push!(indices, par.index - N)
            end
            cur = par
        end

    end

    M
end

function tree2array(tree::Tree,
                    gam::Float64;
                    return_leaf_times::Bool=false)
    nodes = tree.nodes
    _2Nm1 = length(nodes)
    N::Int = (_2Nm1+1) / 2

    times = zeros(N-1)
    leaf_times = zeros(N)

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    for i = reverse(indices)
        cur = tree.nodes[i]
        p = cur.parent
        if i <= N
            self_direction = find(p.children .== cur)[1];
            cur_split = self_direction == 1 ? p.rho : 1-p.rho
            leaf_times[i] = times[p.index-N]*(cur.rhot*cur_split)^gam
            continue
        end

        if p == Nil()
            times[i-N] = (cur.rhot)^gam
        else
            self_direction = find(p.children .== cur)[1];
            cur_split = self_direction == 1 ? p.rho : 1-p.rho
            times[i-N] = times[p.index-N]*(cur.rhot*cur_split)^gam
        end
    end

    I = sortperm(times)
    II = N.+I

    #maps from new to original indices
    sorted_inds = [[1:N], II]

    #maps from original to new indices
    ind_map = zeros(length(sorted_inds))
    ind_map[sorted_inds] = 1:length(sorted_inds)


    Z = zeros(N-1,4)

    for k = 1:length(II)
        i = I[k]
        ii = II[k]
        assert(k+N == ind_map[ii])

        Z[k,3] = times[i]
        l = ind_map[nodes[ii].children[2].index]
        r = ind_map[nodes[ii].children[1].index]
        Z[k,1] = l
        Z[k,2] = r

        l_leaves = l <= N ? 1 : Z[l-N,4]
        r_leaves = r <= N ? 1 : Z[r-N,4]

        Z[k,4] = l_leaves + r_leaves
    end

    states = [nodes[i].state for i = sorted_inds]
    
    if return_leaf_times
        (Z, leaf_times, states, sorted_inds)
    else
        (Z, states, sorted_inds)
    end
end

function serialize(stream, tree::Tree)
    nodes = tree.nodes
    _2Nm1 = length(nodes)

    states = [nodes[i].state for i = 1:_2Nm1] 
    rhos = [nodes[i].rho for i = 1:_2Nm1]
    rhots = [nodes[i].rhot for i = 1:_2Nm1]
   
    n_leaves = [nodes[i].num_leaves for i = 1:_2Nm1]
    n_ancestors = [nodes[i].num_ancestors for i = 1:_2Nm1]

    parents = [ get_parent_index(nodes[i]) for i = 1:_2Nm1]
    childrens = [ get_children_indices(nodes[i]) for i = 1:_2Nm1] 

   
    all_values = (states, rhos, rhots, parents, childrens, n_leaves, n_ancestors)
    serialize(stream, all_values) 
end

function get_parent_index(node::TreeNode)
    p = node.parent
    return p == Nil() ? 0 : p.index
end

function get_children_indices(node::TreeNode)
    c = node.children
    c1 = c[1] == Nil() ? 0 : c[1].index
    c2 = c[2] == Nil() ? 0 : c[2].index

    return (c1,c2)
end

function deserializeTree(stream)

    (states, rhos, rhots, parents, childrens, n_leaves, n_ancestors) = deserialize(stream)

    tree = Tree(states)
    nodes = tree.nodes

    _2Nm1 = length(states)
    for i = 1:_2Nm1
        nodes[i].rho = rhos[i]
        nodes[i].rhot = rhots[i]

        nodes[i].num_leaves = n_leaves[i]
        nodes[i].num_ancestors = n_ancestors[i]

        p_index = parents[i]
        c1_index, c2_index = childrens[i]

        nodes[i].parent = p_index == 0 ? Nil() : nodes[p_index]
        nodes[i].children[1] = c1_index == 0 ? Nil() : nodes[c1_index]
        nodes[i].children[2] = c2_index == 0 ? Nil() : nodes[c2_index]

    end

    tree
end



function total_branch_length(tree::Tree,
                             gam::Float64)
    nodes = tree.nodes
    _2Nm1 = length(nodes)
    N::Int = (_2Nm1+1) / 2

    times = zeros(N-1)

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    total_length = 0.0

    for i = reverse(indices)
        cur = tree.nodes[i]
        p = cur.parent
        if i <= N
            continue
        end

        if p == Nil()
            times[i-N] = (cur.rhot)^gam
            total_length += 1-times[i-N]
        else
            self_direction = find(p.children .== cur)[1];
            cur_split = self_direction == 1 ? p.rho : 1-p.rho
            times[i-N] = times[p.index-N]*(cur.rhot*cur_split)^gam
            total_length += times[p.index-N] - times[i-N]
        end
    end

    return total_length
end

function GenerateUniformTree(N::Int64)

    nodes = Array(TreeNode{Int64}, 2N-1)  
    tree = Tree(nodes)
    for i = 1:N
        cur = TreeNode(i,Nil(),Array(Node,2),i,1,0,1.0,0.5)
        cur.children[1] = Nil()
        cur.children[2] = Nil()
        nodes[i] = cur
        if i > 1
            parent = TreeNode(i+N-1, Nil(),Array(Node,2),i+N-1,1,0,1.0,0.5)
            parent.children[2] = Nil()
            cur.parent = parent
            parent.children[1] = cur
            nodes[i+N-1] = parent

            j = rand() < (i-1)/(2i-3) ? rand(1:i-1) : N+rand(1:i-2)
            InsertIndexIntoTree!(tree, i, j)
        end 
    end
   
    tree 
end
