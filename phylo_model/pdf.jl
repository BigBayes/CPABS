require("phylo_model/phylo_model.jl")
require("tree.jl")
require("utils/probability_util.jl")

using NumericExtensions
import NumericExtensions.evaluate, NumericExtensions.result_type

#require("profile.jl")
#Global pdf
function full_pdf(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState)
    prior(model,model_spec)+likelihood(model,model_spec,data)
end

function prior(model::ModelState,
               model_spec::ModelSpecification)

    Z = model.Z 

    gam = model.gamma
    lambda = model.lambda
    alpha = model.alpha

    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    psi_term = poisson_logpdf(N-1, lambda)

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    mu = ones(2N - 1)
    t = ones(2N - 1)
    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            cur_mu_prop *= cur.rhot
            mu[i] = mu[parent.index]*cur_mu_prop
            t[i] = mu[i]^gam
        else
            mu[i] = cur.rhot^gam
            t[i] = cur.rhot^gam
        end
    end

    #t[1:N] = 0.0

    Tau = zeros(N-1)

    for i = 1:length(tree.nodes)
        cur = tree.nodes[i]
        parent = cur.parent
        if i > N
            left_child = tree.nodes[i].children[2]
            right_child = tree.nodes[i].children[1]
            left_prob = 1.0-tree.nodes[i].rho
            right_prob = tree.nodes[i].rho
    
            N_i = cur.num_leaves
            psi_term += (left_child.num_leaves-1) * log(left_prob)
            psi_term += (right_child.num_leaves-1) * log(right_prob)
            psi_term += -log(N_i) - log(N_i-1)

            N_i = cur.num_leaves
            p_s = 1-2/(N_i+1)

            rhot = tree.nodes[i].rhot
            if rhot == 1.0
                psi_term += log(p_s)
            else
                psi_term += log(1-p_s) + log(p_s) + (p_s-1)*log(rhot)
            end
    
            tau_node = tau(cur)
            tau_t = tau_node == Nil() ? 1.0 : t[tau_node.index]
            Tau[i-N] = tau_t - t[i]

        end
        parent_t = 1.0
        if i != root.index
            parent_t = t[parent.index]
        end

        
    end

    U = zeros(N-1)
    for i=1:length(Z)
        U[Z[i]-N] += 1
    end 
    U = U .- 1
    if any(U .< 0) # There must be at least one point assigned to each cluster
        return -Inf
    end

    Tau = Tau ./ sum(Tau)

    assignment_term = 0.0
    for i = 1:length(U)
        assignment_term += U[i]*log(Tau[i])
    end

    psi_term + assignment_term
end


function likelihood(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

    AA = data.reference_counts
    DD = data.total_counts

    mu_r = data.mu_r
    mu_v = data.mu_v

    eta_term = 0.0
    total_prob = 0.0

    tree = model.tree
    gam = model.gamma
    alpha = model.alpha
    Z = model.Z

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)


    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    num_mutations, num_samples = size(AA)

    t = ones(2N - 1)

    B = zeros(2N - 1, num_samples)

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            t[i] = t[parent.index]*(cur.rhot*cur_mu_prop)^gam

            for s = 1:num_samples
                # the eta variable held by a node is the eta for the right child
                eta_self = self_direction == 1 ? parent.state[s] : (1 - parent.state[s])
                B[i,s] = B[parent.index,s]*eta_self
            end
        else
            t[i] = cur.rhot^gam
            B[i,:] = 1.0
        end

        if i > N
            left_prob = 1.0-tree.nodes[i].rho
            right_prob = tree.nodes[i].rho
            for s = 1:length(tree.nodes[i].state)
                eta_term += logpdf(Beta(alpha*right_prob+1,alpha*left_prob+1), tree.nodes[i].state[s]) 
            end
        end

    end

    for j = 1:num_mutations
        z = Z[j]
        cur_eta = tree.nodes[z].state
        for s = 1:num_samples
            phi_z = B[z,s]*cur_eta[s]
            reference_allele_prob = (1-phi_z)*mu_r[j] + phi_z*mu_v[j]
#            if j == 1
#                println("j=1 s=$s k=$z, phi, ref prob: $phi_z $reference_allele_prob")
#            end
            total_prob += logpdf(Binomial(DD[j,s], reference_allele_prob), AA[j,s])
        end
        
    end

    total_prob + eta_term
end

# Local pdfs for sampling

###################################
###### pdfs for psi updates #######
###################################

#@profile begin
#pdf for updating psi, the tree structure, assumes tree is already pruned at prune_index
function psi_infsites_logpdf(model::ModelState,
                             data::DataState,
                             pruned_index::Int64,
                             path::Array{Int64,1})
    tree = model.tree
    gam = model.gamma
    lambda = model.lambda
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    Z = model.Z
    U = zeros(N-1)
    for i=1:length(Z)
        U[Z[i]-N] += 1
    end 
    U = U .- 1


    #psi_probs = prior_path(tree, pruned_index, path_leaf_index)
    psi_probs = prior_tree(tree, pruned_index)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    total_probs = Float64[]
    tree_states = Array(Int64,0)

    i = 1
    while in(i, subtree_indices)
        i += 1
    end
    root = FindRoot(tree, i)

    indices = GetLeafToRootOrdering(tree, root.index)
    pruned_indices = GetLeafToRootOrdering(tree, pruned_index)

    mu_1 = ones(2N - 1)
    mu_2 = ones(2N - 1)
    mu_3 = ones(2N - 1)

    pruned_node = tree.nodes[pruned_index]
    pruned_parent = pruned_node.parent

    self_direction = find(pruned_parent.children .== pruned_node)[1]

    # We should only be pruning right children
    @assert self_direction == 1

    pruned_mu_prop = pruned_parent.rhot*(1-pruned_parent.rho) 

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    
        num_leaves = cur.num_leaves

        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
            mu_1[i] = mu_1[parent.index]*cur_mu_prop*cur.rhot
            mu_2[i] = mu_1[parent.index]*cur_mu_prop 
        else
            mu_1[i] = cur.rhot
            mu_2[i] = 1.0
        end


    end
    for j = reverse(pruned_indices)
        cur = tree.nodes[j]
        parent = cur.parent    
        if j == pruned_index
            mu_3[parent.index] = parent.rhot
        end

        self_direction = find(parent.children .== cur)[1]
        cur_mu_prop = self_direction == 1 ? parent.rho : 1-parent.rho
        cur_mu_prop *= cur.rhot
        mu_3[j] = mu_3[parent.index]*cur_mu_prop
    end

    t_1 = mu_1.^gam
    t_2 = mu_2.^gam
    t_3 = mu_3.^gam

#    t_1[1:N] = 0.0
#    t_2[1:N] = 0.0
#    t_3[1:N] = 0.0

    # probability of the pruned tree, neglecting the 1/V of the stick lengths for now
    orignal_prob = 0.0


    # cumulative sums of U from leaves to each node
    U_sub = zeros(2N-1)
    for i = indices
        cur = tree.nodes[i]
        if i > N
            l = cur.children[2].index
            r = cur.children[1].index
        
            U_sub[i] = U_sub[l] + U_sub[r] + U[i-N] 
        end
    end

    # Tau, the first ancestor reached from the right child
    Tau = zeros(2N-1)
    for i = reverse(indices)
        cur = tree.nodes[i]
        if i == root.index
            Tau[i] = 0
        end

        if i > N
            l = cur.children[2].index
            r = cur.children[1].index

            Tau[l] = Tau[i]
            Tau[r] = i 
        end 
    end
    for j = reverse(pruned_indices)
        cur = tree.nodes[j]
        if j == pruned_parent.index
            Tau[j] = 0
        end
        if j > N
            l = cur.children[2].index
            r = cur.children[1].index

            Tau[l] = Tau[j]
            Tau[r] = j
        end
    end


    # Distances from each node to its respective tau
    L = zeros(2N-1)
    L_pruned = 0.0

    for i = indices
        if i > N
            tau_index = Tau[i]
            tau_t = tau_index == 0 ? 1.0 : t_1[tau_index]

            L[i] = tau_t - t_1[i]
        end
    end

    for j = pruned_indices
        if j > N
            tau_index = Tau[j]
            tau_t = tau_index == 0 ? 1.0 : t_3[tau_index]

            L[j] = tau_t - t_3[j]
            L_pruned += L[j]
        end
    end

    # cumulative sums of L 
    L_sub = zeros(2N-1)
    for i = indices
        if i > N
            cur = tree.nodes[i]
            l = cur.children[2].index
            r = cur.children[1].index
            L_sub[i] = L[i] + L_sub[l] + L_sub[r]
        end
    end

    sum_L = sum(L[indices])
    sum_U = sum(U)

    for i = indices
        if !( i in path)
            continue
        end
        cur = tree.nodes[i]
        parent = cur.parent 

        pruned_subtree_mu_prop = mu_2[i]
        grafted_subtree_mu_prop = pruned_mu_prop


        subtree_probs = 0.0
        for j = pruned_indices
            if j > N
                L_new = L[j] * pruned_subtree_mu_prop^gam
                subtree_probs += U[j-N]*log(L_new)
            end
        end

        L_new_total = 0.0
        child_indices = GetLeafToRootOrdering(tree, i)
        for j = child_indices
            if j > N 
                tau_index = Tau[j]
                if tau_index == 0
                    tau_t = 1.0
                else
                    tau_t = t_1[tau_index]
                    if tau_t <= t_1[i]
                        tau_t *= grafted_subtree_mu_prop^gam
                    end
                end

                L_new = tau_t - grafted_subtree_mu_prop^gam * t_1[j]
                L_new_total += L_new
                # need to compute the difference here as these terms aren't the same for each i
                subtree_probs += U[j-N]*(log(L_new) - log(L[j]))
            end
        end 

        # one term remains for the "pruned parent"
        # tau of the pruned parent is the same as that of the graftpoint i
        tau_index = Tau[i]
        tau_t = tau_index == 0 ? 1.0 : t_1[tau_index]

        if parent != Nil()
            i_direction = find(parent.children .== cur)[1]
            nu_i = i_direction == 1 ? parent.rho : 1-parent.rho
            pruned_parent_t = (pruned_parent.rhot * nu_i * mu_1[parent.index])^gam
        else
            pruned_parent_t = pruned_parent.rhot^gam
        end

        L_pruned_parent = tau_t - pruned_parent_t

        pruned_parent_term = L_pruned_parent == 0.0 ? -Inf : U[pruned_parent.index-N]*log(L_pruned_parent)
        subtree_probs += pruned_parent_term 

        L_graft = sum_L - L_sub[i] + L_new_total + L_pruned_parent + pruned_subtree_mu_prop^gam * L_pruned
        subtree_probs -= sum_U * log(L_graft)

        push!(total_probs, psi_probs[i] + subtree_probs)
        if isnan(total_probs[end])
            println("psi prob ", psi_probs[i])
            println("subtree_probs ", subtree_probs)
        end

        # node above which to attach
        push!(tree_states, i)                
    end

    (total_probs, tree_states)
end

# Likelihood for mcmc moves for psi
# likelihood for updating psi, the tree structure, assumes tree is already pruned at prune_index
function psi_observation_logpdf(model::ModelState,
                                model_spec::ModelSpecification,
                                data::DataState,
                                pruned_index::Int64,
                                path::Array{Int64,1})

    AA = data.reference_counts
    DD = data.total_counts

    mu_r = data.mu_r
    mu_v = data.mu_v

    tree = model.tree
    gam = model.gamma
    lambda = model.lambda
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    num_mutations, num_samples = size(AA)

    Z = model.Z

    C = [Int64[] for x = 1:maximum(Z)]

    for i = 1:length(Z)
        push!(C[Z[i]], i)
    end
    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    total_probs = Float64[]
    tree_states = Array(Int64,0)

    i = 1
    while in(i, subtree_indices)
        i += 1
    end
    root = FindRoot(tree, i)

    indices = GetLeafToRootOrdering(tree, root.index)
    pruned_indices = GetLeafToRootOrdering(tree, pruned_index)

    # cumulative product of eta variables from root to each node
    Eta = ones(2N-1, num_samples)
    for i = reverse(indices)
        cur = tree.nodes[i]
        if i > N
            l = cur.children[2].index
            r = cur.children[1].index
            for s = 1:num_samples
                eta = cur.state[s]
                Eta[l,s] = (1-eta)*Eta[i,s]
                Eta[r,s] = eta*Eta[i,s]
            end 
        end
    end

    pruned_Eta = ones(2N-1, num_samples)
    for j = reverse(pruned_indices)
        cur = tree.nodes[j]
        if j == pruned_index
            pruned_Eta[j,:] = cur.parent.state
        end
        if j > N
            l = cur.children[2].index
            r = cur.children[1].index
          
            for s = 1:num_samples 
                eta = cur.state[s]
                pruned_Eta[l,s] = (1-eta)*Eta[j,s]
                pruned_Eta[r,s] = eta*Eta[j,s]
            end
        end
    end

    for i = indices
        if !( i in path)
            continue
        end
        cur = tree.nodes[i]

        parent = cur.parent


        subtree_probs = 0.0
        for j = pruned_indices
            if j > N
                cur = tree.nodes[j]
                r = cur.children[1]

                for s = 1:num_samples
                    phi_j = Eta[i,s] * pruned_Eta[j,s] * cur.state[s]

                    for c = C[j]
                        p_j = (1-phi_j)*mu_r[c] + phi_j * mu_v[c]
                        subtree_probs += logpdf(Binomial(DD[c,s], p_j), AA[c,s]) 
                    end
                end
            end
        end


        child_indices = GetLeafToRootOrdering(tree, i)
        for j = child_indices
            if j > N
                cur = tree.nodes[j]
                r = cur.children[1]

                for s = 1:num_samples
                    Eta_l = 1.0 - pruned_Eta[pruned_index,s]
                    phi_j = Eta[j,s] * Eta_l * cur.state[s]

                    old_phi_j = Eta[j,s] * cur.state[s] 
                    
                    for c = C[j]
                        p_j = (1-phi_j)*mu_r[c] + phi_j * mu_v[c]
                        subtree_probs += logpdf(Binomial(DD[c,s], p_j), AA[c,s])
                        # need to subtract "old" terms as terms computed depends on i
                        # (in the pruned subtree, we always compute a term for each node in the pruned tree) 
                        old_p_j = (1 - old_phi_j)*mu_r[c] + old_phi_j * mu_v[c]
                        subtree_probs -= logpdf(Binomial(DD[c,s], old_p_j), AA[c,s])
                    end
                end 
            end
        end 

        # Need a term for the pruned parent
        for s = 1:num_samples
            pruned_parent_index = tree.nodes[pruned_index].parent.index
            phi_j = Eta[i,s] * tree.nodes[pruned_parent_index].state[s]
            for c = C[pruned_parent_index]
                p_j = (1-phi_j)*mu_r[c] + phi_j * mu_v[c]
                subtree_probs += logpdf(Binomial(DD[c,s], p_j), AA[c,s])
            end

        end


        push!(total_probs, subtree_probs)
        if isnan(total_probs[end])
            println("subtree_probs ", subtree_probs)
        end
        push!(tree_states, i)                
    end

    (total_probs, tree_states)
end


function prior_tree(tree::Tree,
                    pruned_index::Int)
    N::Int = (length(tree.nodes) + 1) / 2

    prob_diffs = -Inf*ones(2N - 1)
    prob_left = -Inf*ones(2N - 1)
    prob_right = -Inf*ones(2N - 1)
    rhot_prob = -Inf*ones(2N-1)

    subtree_indices = GetSubtreeIndicies(tree, pruned_index)

    i = 1
    while in(i,subtree_indices)
        i += 1
    end
    root = FindRoot(tree, i)

    # indices will contain a leaves to root ordering of nodes in the pruned tree
    indices = GetLeafToRootOrdering(tree, root.index)

    pruned_indices = GetLeafToRootOrdering(tree, pruned_index)
    pruned_node = tree.nodes[pruned_index]
    pruned_parent = pruned_node.parent

    self_direction = find(pruned_parent.children .== pruned_node)[1]

    pruned_mu_prop = self_direction == 1 ? 1-pruned_parent.rho : (pruned_parent.rho) 
    pruned_leaves = pruned_node.num_leaves
    
    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    
        num_leaves = cur.num_leaves

        if i == root.index
            parent_rhot_prob = 0.0
            cur_prob = 0.0
        else
            parent_rhot_prob = rhot_prob[parent.index]
            self_direction = find(parent.children .== cur)[1]
            cur_prob = self_direction == 1 ? prob_right[parent.index] : prob_left[parent.index]
        end

        if i > N
            N_l = cur.children[2].num_leaves
            N_r = cur.children[1].num_leaves
            N_i = cur.num_leaves
            assert(N_l + N_r == N_i)

            ps_new = 1-2/(N_i+pruned_leaves+1)
            ps_old = 1-2/(N_i+1)

            rhot_prob[i] = cur.rhot == 1.0 ? log(ps_new)-log(ps_old) : 
                                (log(1-ps_new) + log(ps_new) + (ps_new-1)*log(cur.rhot) -
                                 log(1-ps_old) - log(ps_old) - (ps_old-1)*log(cur.rhot))
            rhot_prob[i] += parent_rhot_prob

            prob_right[i] = pruned_leaves * log(cur.rho) -
                            log(N_i+pruned_leaves) - log(N_i+pruned_leaves-1) +
                            log(N_i) + log(N_i-1) +
                            cur_prob
            prob_left[i] = pruned_leaves * log(1-cur.rho) - 
                           log(N_i+pruned_leaves) - log(N_i+pruned_leaves-1) +
                           log(N_i) + log(N_i-1) +
                           cur_prob
        end

        total_leaves = num_leaves + pruned_leaves
        p_s = 1-2/(total_leaves + 1)
        rhot = pruned_parent.rhot
        cur_rhot_prob = rhot == 1.0 ? log(p_s) : log(1-p_s) + log(p_s) + (p_s-1)*log(rhot) 

        prob_diffs[i] = cur_prob + parent_rhot_prob + cur_rhot_prob +
                        (num_leaves - 1) * log(pruned_mu_prop) -
                        log(total_leaves) - log(total_leaves-1) 
#        prob_above[i] = num_leaves * log(pruned_mu_prop)
#        prob_below[i] = pruned_leaves * log(mu_1[i])

#        if i != root.index
#            cum_prob_above[i] = cum_prob_above[parent.index] + prob_above[i]
#        end
    end


    prob_diffs
end

function prune_graft_logprobs(input_model::ModelState,
                              model_spec::ModelSpecification,
                              data::DataState,
                              prune_index::Int64)

    model = copy(input_model)
    tree = model.tree
    PruneIndexFromTree!(tree, prune_index)
    subtree_indices = GetSubtreeIndicies(model.tree, prune_index)
    i = 1
    while in(i,subtree_indices)
        i += 1
    end
    root = FindRoot(tree, i)
    path = GetLeafToRootOrdering(tree, root.index)

    priors = Float64[]
    likelihoods = Float64[]

    for graft_index = path
        InsertIndexIntoTree!(tree, prune_index, graft_index)
        prior_term = prior(model, model_spec)
        likelihood_term = likelihood(model, model_spec, data)

        push!(priors, prior_term)
        push!(likelihoods, likelihood_term)
        PruneIndexFromTree!(tree, prune_index)
    end

    return (priors, likelihoods)
end

###################################
###### pdfs for nu, nutd updates #
###################################

function nu_tilde_splits(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r, xi; node="r")

    l1 = node == "r" ? nutd_l : 1.0
    r1 = node == "r" ? 1.0 : nutd_r

    p_z1 = p_z_given_nu_nutd(nu, l1, r1, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r)
    p_z2 = p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r)
 
    f1 = log(xi) + p_z1

    f2 = node == "r" ? log(1-xi)+log(xi)+(xi-1)*log(nutd_r) + p_z2 :
                       log(1-xi)+log(xi)+(xi-1)*log(nutd_l) + p_z2

#    println("xi: $xi")
#    println("nu: $nu")
#    println("nutd_l: $nutd_l")
#    println("nutd_r: $nutd_r")
#    println("splits: $f1 $f2")
#
#    if isnan(f1)
#        println("f1")
#    elseif isnan(f2)
#        println("p_z2: $p_z2")
#        if node == "r"
#            println("r.  xi: $xi, nutd_r: $nutd_r")
#        else
#            println("l.  xi: $xi, nutd_l: $nutd_l")
#        end
#    end

    return [f1, f2] 
end

function root_nu_tilde_splits(nutd, gam, U, U_i, u, K_i, C_i, D, Pki_i, xi)

    p_z1 = p_z_given_nutd(1.0, gam, U, U_i, u, K_i, C_i, D, Pki_i)
    p_z2 = p_z_given_nutd(nutd, gam, U, U_i, u, K_i, C_i, D, Pki_i)
 
    f1 = log(xi) + p_z1

    f2 = log(1-xi)+log(xi)+(xi-1)*log(nutd) + p_z2

    return [f1, f2] 
end

function nu_logpdf(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r,
                    eta, alpha, N_l, N_r)

    p_z = p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r)
    p_eta = p_eta_given_nu(eta, nu, alpha)
    p_nu = p_nu_Nl_Nr(nu, N_l, N_r)

    p_z + p_eta + p_nu
end

function p_nu_Nl_Nr(nu, N_l, N_r)
    (N_r-1)log(nu) + (N_l-1)log(1-nu) - log(N_l+N_r) - log(N_l+N_r-1)
end

function p_eta_given_nu(eta, nu, alpha)
    sum(logpdf(Beta(alpha*nu+1, alpha*(1-nu)+1), eta))
end

function p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, Pki_l, Pki_r)

    if nu == 0.0 || nu == 1.0 || nutd_r == 0.0 || nutd_l == 0.0
        return -Inf
    end

    log_pdf = gam*U_r*log(nu) + gam*U_l*log(1-nu) + gam*U_l*log(nutd_l)+gam*U_r*log(nutd_r) - 
              U*log( (nu*nutd_r)^gam * C_r + ((1-nu)*nutd_l)^gam * C_l + D )

    rrg = ((1-nu)*nutd_l)^gam

    for c = Pki_l
        if c != Nil()
            i = c.index
            if 1-rrg*K_l[i] == 0.0
                log_pdf += -Inf
            else 
                log_pdf += u[i]*log(1-rrg*K_l[i])
            end
        end
    end

    rrg = (nu*nutd_r)^gam
    for c = Pki_r
        if c != Nil()
            i = c.index 
            if 1-rrg*K_r[i] == 0.0
                log_pdf += -Inf
            else 
                log_pdf += u[i]*log(1-rrg*K_r[i])
            end
        end
    end
    
    log_pdf
end

function p_z_given_nutd(nutd, gam, U, U_i, u, K_i, C_i, D, Pki_i)

    if nutd == 0.0 
        return -Inf
    end

    log_pdf = gam*U_i*log(nutd) - U*log(nutd^gam * C_i + D)


    rrg = nutd^gam
    for c = Pki_i
        if c != Nil()
            i = c.index
            if 1-rrg*K_i[i] == 0.0
                log_pdf += -Inf
            else 
                log_pdf += u[i]*log(1-rrg*K_i[i])
            end
        end
    end

    log_pdf
end

###################################
###### pdfs for eta updates #######
###################################

function eta_logpdf(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState;
                    Temp::Float64 = 1.0)

    tree = model.tree
    N::Int = (length(tree.nodes) + 1) / 2

    alpha = model.alpha
    Z = model.Z

    C = [Int64[] for x = 1:maximum(Z)]
    for i = 1:length(Z)
        push!(C[Z[i]], i)
    end

    AA = data.reference_counts
    DD = data.total_counts
    mu_r = data.mu_r
    mu_v = data.mu_v

    M, S = size(AA)

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    log_pdf = 0.0

    Phi = zeros(2N - 1, S)

    for j = reverse(indices)



        if j > N
            cur = tree.nodes[j]
            eta = cur.state
            nu_r = cur.rho

            log_pdf += (alpha*nu_r)*sum(log(eta)) + (alpha*(1-nu_r))*sum(log(1 .- eta))

 
            phi = ones(size(eta))
            An = GetAncestors(tree, cur.index)
            for a = An
                p = a.parent
                if p != Nil()
                    a_direction = find(p.children .== a)[1]
                    eta_a = a_direction == 1 ? p.state : 1.0 .- p.state
                    phi .*= eta_a
                end
            end
            phi .*= eta
            for i = C[j]
                r = (1.-phi)*mu_r[i] + phi*mu_v[i]
                logr = log(r)
                log1mr = log(1.-r) 
                for s = 1:length(eta)
                    log_pdf += (AA[i,s]*log(r[s]) + (DD[i,s]-AA[i,s])*log(1-r[s]))/Temp
                end
            end

        end 
    end

    log_pdf
end

function eta_log_gradient(model::ModelState,
                          model_spec::ModelSpecification,
                          data::DataState;
                          Temp::Float64 = Temp)
    tree = model.tree
    N::Int = (length(tree.nodes) + 1) / 2

    alpha = model.alpha
    Z = model.Z

    C = [Int64[] for x=1:maximum(Z)]

    for i = 1:length(Z)
        push!(C[Z[i]], i)
    end

    AA = data.reference_counts
    DD = data.total_counts
    mu_r = data.mu_r
    mu_v = data.mu_v

    (M,S) = size(AA)

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    gradient = zeros(N-1, S)
    phi = ones(N-1, S)
    for j = indices
        if j > N
            cur = tree.nodes[j]
          
            An = GetAncestors(tree, cur.index)
            for a = An
                p = a.parent
                if p != Nil()
                    a_direction = find(p.children .== a)[1]
                    eta_a = a_direction == 1 ? p.state : 1.0 .- p.state
                    phi[j-N, :] .*= eta_a'

                    # return gradient = 0 if prior is 0
                    if any(eta_a .== 0.0) || any(eta_a .== 1.0)
                        return gradient
                    end
                end
            end
            phi[j-N,:] .*= cur.state'
        end
    end

    for k = indices
        if k > N
            
            cur = tree.nodes[k]
            parent = cur.parent

            nu_k = parent == Nil() ? 1.0 : cur.parent.rho
            eta_k = cur.state 

            for s = 1:S
                gradient[k-N,s] += (alpha*nu_k)/eta_k[s] - (alpha*(1-nu_k))/(1-eta_k[s])
            end         

            # right subtree 
            De = GetDescendants(tree, cur.children[1].index)
            #println("De($(cur.children[1].index)): $([x.index for x in De])")
            for d = De
                j = d.index
                if j > N
                    phi_j_k = phi[j-N,:]
                    for s = 1:S
                            
                        phi_j_k[s] = phi_j_k[s] == 0.0 ? 0.0 : phi_j_k[s] / eta_k[s]
                    end
                    for i = C[j]
                        r = (1 .- phi[j-N,:]).*mu_r[i] + phi[j-N,:].*mu_v[i]
                        for s = 1:S
                            gradient[k-N,s] += (AA[i,s]*(mu_v[i]-mu_r[i])*phi_j_k[s]/r[s] +
                                               (DD[i,s] - AA[i,s])*(mu_r[i]-mu_v[i])*phi_j_k[s]/(1-r[s]))/Temp
                        end
                    end
                end
            end
            # left subtree 
            De = GetDescendants(tree, cur.children[2].index)
            for d = De
                j = d.index
                if j > N
                    phi_j_k = phi[j-N,:]
                    for s = 1:S
                        phi_j_k[s] = phi_j_k[s] == 0.0 ? 0.0 : phi_j_k[s] / (1 - eta_k[s])
                    end
                    for i = C[j]
                        r = (1 .- phi[j-N,:]).*mu_r[i] + phi[j-N,:].*mu_v[i]
                        for s = 1:S
                            gradient[k-N,s] -= (AA[i,s]*(mu_v[i]-mu_r[i])*phi_j_k[s]/r[s] +
                                               (DD[i,s] - AA[i,s])*(mu_r[i]-mu_v[i])*phi_j_k[s]/(1-r[s]))/Temp
                        end
                    end
                end
            end

            # current node
            phi_k_k = phi[k-N,:]
            for s = 1:S
                phi_k_k[s] = phi_k_k[s] == 0.0 ? 0.0 : phi_k_k[s] / eta_k[s]
            end
            for i = C[k]
                r = (1 .- phi[k-N,:]).*mu_r[i] + phi[k-N,:].*mu_v[i]
                for s = 1:S
                    gradient[k-N,s] += (AA[i,s]*(mu_v[i]-mu_r[i])*phi_k_k[s]/r[s] +
                                       (DD[i,s] - AA[i,s])*(mu_r[i]-mu_v[i])*phi_k_k[s]/(1-r[s]))/Temp
                end
            end
            

        end 

    end
    gradient'[:]
end

###################################
##### pdf for Z assignment update #
###################################

function z_logpdf(model::ModelState,
                  model_spec::ModelSpecification,
                  data::DataState,
                  index::Int64,
                  U::Vector{Int64},
                  times::Vector{Float64},
                  Tau::Vector{Int64},
                  phi::Matrix{Float64};
                  Temp::Float64 = 1.0)

    tree = model.tree
    N::Int = (length(tree.nodes) + 1) / 2

    alpha = model.alpha
    Z = model.Z

    current_k = Z[index]

    # We cannot move this mutation as it would create an empty cluster
    if U[current_k-N] == 0
        result = -Inf*ones(N-1)
        result[current_k-N] = 0.0
        return result
    end

    AA = data.reference_counts
    DD = data.total_counts
    mu_r = data.mu_r
    mu_v = data.mu_v

    (M,S) = size(AA)

    log_pdf = zeros(N-1)

    for k = N+1:2N-1
        tau_t = Tau[k] == 0 ? 1.0 : times[Tau[k]]
        v_k = tau_t - times[k]

        n_k = U[k-N]+1
        q_k = (1 + (n_k-1)*v_k)/n_k
        log_pdf[k-N] = log(q_k)

        for s = 1:S
            r_ks = (1-phi[k,s])*mu_r[index] + phi[k,s]*mu_v[index]
#            if index == 1
#                println("j=1 s=$s k=$k, phi, ref prob (fast): $(phi[k,s]) $r_ks")
#            end
            log_pdf[k-N] += (AA[index,s]*log(r_ks) + (DD[index,s] - AA[index,s])*log(1-r_ks))/Temp
        end

    end
    log_pdf 
end

###################################
###### pdfs for vardim updates ####
###################################
#@profile begin


function vardim_local_logpdf(model::ModelState,
                             data::DataState,
                             relevant_pairs::Array{Int64, 2},
                             latent_effects::Array{Array{Float64, 2}, 1},
                             observed_effects::Array{Float64, 2},
                             mixture_component_index::Int64,
                             u::Int64,
                             node_index::Int64,
                             effective_lambda::Float64,
                             w_old::Float64,
                             w_new::Float64,
                             w_is_auxiliary::Bool)

    (K,K) = size(model.weights)
    sigma = model.w_sigma
    if u == model.tree.nodes[node_index].state - 1
        if w_is_auxiliary
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            # W_local_logpdf includes the prior term for w_old as if it were normal,
            # need to subtract it out here
            logprob += aug_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    elseif u == model.tree.nodes[node_index].state
        if w_is_auxiliary
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_old)

            logprob += aug_logpdf( w_new, model.nu) - normal_logpdf(w_old, sigma)
        else
            logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                      observed_effects, w_old, w_new)
        end
    else
        logprob = W_local_logpdf(model,data,relevant_pairs, latent_effects[mixture_component_index],
                  observed_effects, w_old, w_new)
    end

    logprob + poisson_logpdf(u, effective_lambda)
end

function vardim_local_splits(model::ModelState,
                             model_spec::ModelSpecification,
                             data::DataState,
                             relevant_pairs::Array{Int64, 2},
                             latent_effects::Array{Array{Float64, 2}, 1},
                             observed_effects::Array{Float64, 2},
                             u::Array{Int64, 1},
                             weight_index_pointers::Array{Int64, 1},
                             node_index::Int64,
                             effective_lambda::Float64,
                             constant_terms::Array{Float64,1},
                             w_old::Float64,
                             w_new::Float64,
                             w_is_auxiliary::Array{Bool, 1})

    num_components = length(latent_effects)
    logprobs = zeros(num_components)

    for mixture_component = 1:num_components
        logprob = vardim_local_logpdf(model, data, relevant_pairs, latent_effects,
                                      observed_effects, mixture_component, 
                                      u[mixture_component], node_index,
                                      effective_lambda, w_old, w_new,
                                      w_is_auxiliary[mixture_component]) 
        logprobs[mixture_component] = logprob

    end
    logprobs += vardim_multiplier_terms(model_spec, u)
    logprobs + constant_terms
end

function vardim_local_sum(model::ModelState,
                          model_spec::ModelSpecification,
                          data::DataState,
                          relevant_pairs::Array{Int64, 2},
                          latent_effects::Array{Array{Float64, 2}, 1},
                          observed_effects::Array{Float64, 2},
                          u::Array{Int64, 1},
                          weight_index_pointers::Array{Int64, 1},
                          node_index::Int64,
                          effective_lambda::Float64,
                          constant_terms::Array{Float64,1},
                          w_old::Float64,
                          w_new::Float64,
                          w_is_auxiliary::Array{Bool, 1})
    logprobs = vardim_local_splits(model, model_spec, data, relevant_pairs, latent_effects,
                                   observed_effects, u, weight_index_pointers,
                                   node_index, effective_lambda, constant_terms, 
                                   w_old, w_new, w_is_auxiliary)
    logsumexp(logprobs)
end



function vardim_logpdf(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       observed_effects::Array{Float64, 2},
                       u::Int64,
                       L::Int64,
                       aug_k::Int64,
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64)


    tmodel = copy(model)
    tmodel.tree.nodes[node_index].state = L+1
    Z = ConstructZ(tmodel.tree)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + L - 1
    new_index = end_index + 1 

    W = get_W(model, model_spec)

    if u < L
        assert( aug_k < new_index)
        W[:,new_index] = 0.0
        W[new_index,:] = 0.0
        W[aug_k,:] = 0.0
        W[:,aug_k] = 0.0
    elseif L == u
        W[:,new_index] = 0.0
        W[new_index,:] = 0.0
    elseif u > L
        assert(aug_k == new_index)
    end

    if model_spec.symmetric_W
        W = copy(W)
        symmetrize!(W)
    end

    YY = data.Ytrain
    (N,N) = size(YY[1])

    logprob = 0.0

    gam = model.gamma
    lambda = model.lambda
    tree = model.tree
    _2Nm1 = length(tree.nodes)
    N::Int = (_2Nm1+1)/2

    ZWZ = Z*W*Z' 
    effects = ZWZ + observed_effects

    logprob = 0.0
    gradient = zeros(size(YY[1]))
    for s = 1:length(YY)
        Y = YY[s]
        LL = broadcast(log_logit, effects, Y)
        LL_grad = broadcast(log_logistic_dx, effects, Y)
        missing_data = find(Y .< 0)
        LL[missing_data] = 0.0
        LL_grad[missing_data] = 0.0
      
        logprob += sum(LL) 
        gradient += LL_grad 
    end

    ZLZ = Z'*gradient*Z

    if model_spec.symmetric_W
        ZLZ_L = tril(ZLZ)
        ZLZ_L[diagind(ZLZ)] = 0.0
        ZLZ += ZLZ_L'
        symmetrize!(ZLZ)
    end

    if model_spec.positive_W
        gradient = ZLZ.*W
    else
        gradient = ZLZ
    end

    W_logpdf = model_spec.W_logpdf
    W_logpdf_gradient = model_spec.W_logpdf_gradient
    (K,K) = size(W)
    sigma = model.w_sigma

    W_logprob = W_logpdf(W, sigma)    

    if model_spec.diagonal_W
        logprob += sum(diag(W_logprob)) 
    elseif model_spec.symmetric_W
        logprob += sum(triu(W_logprob))
    else
        logprob += sum(W_logprob)
    end

    W_grad = W_logpdf_gradient(W, sigma)

    if model_spec.positive_W
        gradient += W_grad.*W
    else
        gradient += W_grad
    end

    logprob += poisson_logpdf(u, effective_lambda)

    return logprob, gradient 
end

function vardim_splits(model::ModelState,
                       model_spec::ModelSpecification,
                       data::DataState,
                       observed_effects::Array{Float64, 2},
                       U::Array{Int64, 1},
                       weight_index_pointers::Array{Int64, 1},
                       node_index::Int64,
                       effective_lambda::Float64)
    logprobs = zeros(length(U))
    gradients = cell(length(U))
    assert(length(U) > 1)
    L = U[end-1] 

    start_index = weight_index_pointers[node_index]
    end_index = start_index + L - 1
    new_index = end_index + 1 

    for u_ind = 1:length(U)
        u = U[u_ind]
          
        if u < L 
            aug_k = start_index + u_ind - 1
        elseif u == L
            aug_k = 0
        elseif u > L
            aug_k = new_index
        end

        logprobs[u_ind], gradients[u_ind] = 
                              vardim_logpdf(model, model_spec, data,
                              observed_effects, u, L, aug_k, 
                              weight_index_pointers, node_index,
                              effective_lambda)

    end
    logprobs += vardim_multiplier_terms(model_spec, U)
    logprobs, gradients
end

function vardim_sum(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState,
                    observed_effects::Array{Float64, 2},
                    U::Array{Int64, 1},
                    weight_index_pointers::Array{Int64, 1},
                    node_index::Int64,
                    effective_lambda::Float64)

    logprobs, gradients = vardim_splits(model, model_spec, data, observed_effects,
                                        U, weight_index_pointers, node_index,
                                        effective_lambda)

    logprob = logsumexp(logprobs)
    gradient = logsumexp_d_dx(logprobs, gradients)

    logprob, gradient
end

function vardim_multiplier_terms(model_spec::ModelSpecification,
                                 u::Array{Int64,1})
    
    if length(u) == 1 # the case L < 0
        return [log(model_spec.rrj_jump_probabilities[1])]
    end

    logprobs = zeros(length(u))
    logprobs += log(model_spec.rrj_jump_probabilities[u .- u[1] .+ 1])
    # for u with K < L == u[end-1], we need the 1/L term
    logprobs -= [u[x] < u[end-1] ? log(u[end-1]) : 0.0 for x = 1:length(u)]

    logprobs
end

#end #profile
# Utility Functions

# Assumes model.augmented_weights is up to date!
function compute_constant_terms(model::ModelState,
                                model_spec::ModelSpecification,
                                data::DataState,
                                relevant_pairs::Array{Array{Int64, 2}, 2},
                                latent_effects::Array{Array{Float64, 2}, 1},
                                observed_effects::Array{Float64, 2},
                                U::Array{Int64, 1},
                                weight_index_pointers::Array{Int64, 1},
                                node_index::Int64,
                                effective_lambda::Float64,
                                total_logprob::Float64)

    assert(length(U) > 0)
    logprobs = total_logprob * ones(length(U))

    current_ind = length(U) - 1

    current_u = U[current_ind]
    W = get_W(model, model_spec) 

    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1


    for u_ind = 1:length(U)
        u = U[u_ind]

        aug_k = u < current_u ? start_index + u_ind - 1 : u > current_u ? new_index : 0

        prob_tuple = adjust_model_logprob(model, model_spec, data, relevant_pairs,
                         latent_effects[current_ind], latent_effects[u_ind],
                         observed_effects, aug_k, current_u, u, weight_index_pointers,
                         node_index, effective_lambda)


        logprobs[u_ind] += sum(prob_tuple)
#        if model_spec.debug
#            println("u_ind, prob_tuple: ", (u_ind, prob_tuple))
#            println("u_ind, u, current_u, aug_k: ", (u_ind, u, current_u, aug_k))
#        end
    end 
    logprobs
end



#end #profile

