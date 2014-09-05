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
            mu[i] = mu[parent.index]*cur_mu_prop
            t[i] = t[parent.index]*(cur.rhot*cur_mu_prop)^gam
        else
            t[i] = cur.rhot^gam
        end
    end

    #t[1:N] = 0.0

    eta_term = 0.0
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

            for s = 1:length(tree.nodes[i].state)
                eta_term += logpdf(Beta(alpha*right_prob,alpha*left_prob), cur.state[s]) 
            end
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
    if any(U .< 0) # There must be one point assigned to each cluster
        return -Inf
    end

    Tau = Tau ./ sum(Tau)

    assignment_term = 0.0
    for i = 1:length(U)
        assignment_term += U[i]*log(Tau[i])
    end

    psi_term + eta_term + assignment_term
end


function likelihood(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

    AA = data.reference_counts
    DD = data.total_counts

    mu_r = data.mu_r
    mu_v = data.mu_v

    total_prob = 0.0
    total_test_prob = 0.0

    tree = model.tree
    gam = model.gamma
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

    end

    for j = 1:num_mutations
        z = Z[j]

        for s = 1:num_samples
            phi_z = B[z,s]
            reference_allele_prob = (1-phi_z)*mu_r[j] + phi_z*mu_v[j]
            total_prob += logpdf(Binomial(DD[j,s], reference_allele_prob), AA[j,s])
        end
    end

    total_prob
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

    pruned_mu_prop = self_direction == 1 ? pruned_parent.rho : (1-pruned_parent.rho) 
    pruned_mu_prop *= pruned_parent.rhot

    for i = reverse(indices)
        cur = tree.nodes[i]
        parent = cur.parent    
        num_leaves = cur.num_leaves

        cur_mu_prop = cur.rhot
        if i != root.index
            self_direction = find(parent.children .== cur)[1]
            cur_mu_prop *= self_direction == 1 ? parent.rho : 1-parent.rho
            mu_1[i] = mu_1[parent.index]*cur_mu_prop
        else
            mu_1[i] = cur.rhot
        end

        mu_2[i] = mu_1[i] * pruned_mu_prop 

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
    for j = reverse(indices)
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
        tau_index = Tau[i]
        tau_t = tau_index == 0 ? 1.0 : t_1[tau_index]

        L[i] = tau_t - t_1[i] 
    end

    for j = pruned_indices
        tau_index = Tau[j]
        tau_t = tau_index == 0 ? 1.0 : t_3[tau_index]

        L[j] = tau_t - t_3[j]
        L_pruned += L[j]
    end

    # cumulative sums of L from leaves to each node 
    L_sub = zeros(2N-1)
    for i = indices
        cur = tree.nodes[i]
        if i > N
            l = cur.children[2].index
            r = cur.children[1].index
            L_sub[i] = L[i] + L_sub[l] + L_sub[r]
        else
            L_sub[i] = L[i]
        end
    end

    sum_L = sum(L)
    sum_U = sum(U)

    for i = indices
        if !( i in path)
            continue
        end
        cur = tree.nodes[i]


        pruned_subtree_mu_prop = mu_1[i]
        grafted_subtree_mu_prop = pruned_mu_prop

        L_graft = sum_L - L_sub[i] + grafted_subtree_mu_prop^gam * L_sub[i] + pruned_subtree_mu_prop^gam * L_pruned

        subtree_probs = 0.0
        for j = pruned_indices
            if j > N
                L_new = L[j] * pruned_subtree_mu_prop^gam
                subtree_probs += U[j-N]*log(L_new)
            end
        end

        child_indices = GetLeafToRootOrdering(tree, i)
        for j = child_indices
            if j > N
                L_new = L[i] + pruned_subtree_mu_prop^gam * (t_1[i] - t_1[j])
                # need to compute the difference here as these terms aren't the same for each i
                subtree_probs += U[j-N]*(log(L_new) - log(L[j]))
            end
        end 

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
        elseif j > N
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

        subtree_probs = 0.0
        for j = pruned_indices
            if j > N
                cur = tree.nodes[j]
                r = cur.children[1]

                for s = 1:num_samples
                    phi_j = Eta[i,s] * pruned_Eta[j,s]

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
                    phi_j = Eta[j,s] * Eta_l 

                    old_phi_j = Eta[j,s] 
                    
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

    pruned_mu_prop = self_direction == 1 ? pruned_parent.rho : (1-pruned_parent.rho) 
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


###################################
###### pdfs for nu, nutd updates #
###################################

function nu_tilde_r_splits(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, P_l, P_r, xi; node="r")

    p_z = p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, P_l, P_r)
 
    f1 = log(xi) + p_z

    f2 = node == "r" ? log(1-xi)+log(xi)+(xi-1)*log(nutd_r) + p_z :
                       log(1-xi)+log(xi)+(xi-1)*log(nutd_l) + p_z
    return [f1, f2] 
end


function nu_logpdf(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, P_l, P_r,
                    eta, alpha, N_l, N_r)

    p_z = p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, P_l, P_r)
    p_eta = p_eta_given_nu(eta, nu, alpha)
    p_nu = p_nu_Nl_Nr(nu, N_l, N_r)

    p_z + p_eta + p_nu
end

function p_nu_Nl_Nr(nu, N_l, N_r)
    (N_r-1)log(nu) + (N_l-1)log(1-nu) - log(N_l+N_r) - log(N_l+N_r-1)
end

function p_eta_given_nu(eta, nu, alpha)
    S = length(eta)

    log_pdf = S*lbeta(alpha*nu, alpha*(1-nu))

    leta_1 = sum(log(eta))
    leta_2 = sum(log(1-eta))

    logpdf += (alpha*nu-1)*leta_1 + (alpha*(1-nu)-1)*leta_2
    logpdf
end

function p_z_given_nu_nutd(nu, nutd_l, nutd_r, gam, U, U_l, U_r, u, K_l, K_r, C_l, C_r, D, P_l, P_r)
    log_pdf = gam*U_r*log(nutd_r*nu) + gam*U_l*log(nutd_l*(1-nu)) + gam*U_l*log(nutd_l)+gam*U_r*log(nutd_r) - 
              U*log( (nu*nutd_r)^gam * C_r + ((1-nu)*nutd_l)^gam * C_l + D )

    rrg = ((1-nu)*nutd_l)^gam
    for i = P_l
        log_pdf += u[i]*log(1-rrg*K_l[i])
    end

    rrg = (nu*nutd_r)^gam
    for i = P_r
        log_pdf += u[i]*log(1-rrg*K_r[i])
    end
    
    log_pdf
end


###################################
###### pdfs for eta updates #######
###################################

function eta_logpdf(model::ModelState,
                    model_spec::ModelSpecification,
                    data::DataState)

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

    root = FindRoot(tree, 1)
    indices = GetLeafToRootOrdering(tree, root.index)

    log_pdf = 0.0

    for j = indices
        if j > N
            cur = tree.nodes[j]
            eta = cur.state
            nu_r = cur.parent == Nil() ? 1.0 : cur.parent.rho

            log_pdf += (alpha*nu_r-1)*sum(log(eta)) + (alpha*(1-nu_r)-1)*sum(log(1 .- eta))
          
            phi = ones(length(eta))
            An = GetAncestors(cur)
            for a = An
                p = a.parent
                if p != Nil()
                    a_direction = find(p.children .== a)[1]
                    eta_a = a_direction == 1 ? a.state : 1.0 .- a.state
                    phi .*= eta_a
                end
            end
            for i = C[j]
                r = (1.-phi)*mu_r[i] + phi*mu_v[i]
                logr = log(r)
                log1mr = log(1.-r) 
                for s = 1:length(eta)
                    log_pdf += AA[i,s]*log(r[s]) + (DD[i,s]-AA[i,s])*log(1-r[s])
                end
            end

        end 

    end

    log_pdf
end

function eta_log_gradient(model::ModelState,
                          model_spec::ModelSpecification,
                          data::DataState)
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
          
            An = GetAncestors(cur)
            for a = An
                p = a.parent
                if p != Nil()
                    a_direction = find(p.children .== a)[1]
                    eta_a = a_direction == 1 ? a.state : 1.0 .- a.state
                    phi[j-N, :] .*= eta_a
                end
            end
        end
    end
   
    for k = indices
        if k > N
            
            cur = tree.nodes[k]
            parent = cur.parent
            nu_k = parent == Nil() ? 1.0 : cur.parent.rho
            self_direction = find(parent.children .== cur)[1]
            eta_p = parent.eta
            eta_k = self_direction == 1 ? eta_p : 1.-eta_p 
    
            for s = 1:S
                gradient[k-N,s] += (alpha*nu_k-1)/eta_k[s] + (alpha*(1-nu_k)-1)/(1-eta_k[s])
            end         

 
            De = GetDescendants(cur)
            for j = De
                p = j.parent
                phi_j_k = phi[j,:] ./ eta_k
                for i = C[j]
                    r = (1 .- phi[k,:]).*mu_r[j] + phi[k,:].*mu_v[j]
                    for s = 1:S
                        gradient[k-N,s] +=  AA[i,s]*(mu_v[i]-mu_r[i])*phi_j_k/r[s] +
                                            (DD[i,s] - AA[i,s])*(mu_r[i]-mu_v[i])*phi_j_k/(1-r[s])
                    end
                end
            end
        end 

    end
    gradient
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
                  Tau::Vector{TreeNode},
                  phi::Vector{Float64})

    tree = model.tree
    N::Int = (length(tree.nodes) + 1) / 2

    alpha = model.alpha
    Z = model.Z

    current_k = Z[index]

    # We cannot move this mutation as it would create an empty cluster
    if U[current_k-N] == 0
        result = -Inf*ones(N-1)
        result[current_k] = 0.0
        return result
    end

    AA = data.reference_counts
    DD = data.total_counts
    mu_r = data.mu_r
    mu_v = data.mu_v

    (M,S) = size(AA)

    An = GetAncestors(cur)

    log_pdf = zeros(N-1)

    for k = N+1:2N-1
        v_k = times[Tau[k].index] - times[k]

        n_k = U[k-N]
        q_k = (1 + (n_k-1)*v_k)/n_k
        log_pdf[k] = log(q_k)

        for s = 1:S
            r_ks = (1-phi[k,s])*mu_r[index] + phi[k,s]*mu_v[index]
            log_pdf[k] += AA[index,s]*log(r_ks) + (DD[index,s] - AA[index,s])*log(1-r_ks) 
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

# Adjust model_logprob when L != u
# Assumes model.augmented_weights is up to date!
function adjust_model_logprob(model::ModelState,
                              model_spec::ModelSpecification,
                              data::DataState,
                              relevant_pairs::Array{Array{Int64, 2}, 2},
                              current_latent_effects::Array{Float64, 2},
                              latent_effects::Array{Float64, 2},
                              observed_effects::Array{Float64, 2},
                              aug_k::Int,
                              u::Int,
                              L::Int,
                              weight_index_pointers::Array{Int64, 1},
                              node_index::Int64,
                              effective_lambda::Float64)

    if L < 0
        return (0.0, 0.0, 0.0)
    end

    YY = data.Ytrain

    augW = model.augmented_weights
    aug_inds, _ = get_augmented_submatrix_indices(augW, node_index, 1)

    W = augW[aug_inds, aug_inds]
    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1

    aug_terms = 0.0
    prior_terms = 0.0
    likelihood_terms = 0.0
    if model_spec.diagonal_W
        if L < u
            removed_ind = aug_k
            aug_terms += aug_logpdf(W[new_index,new_index], model.nu)
            aug_terms += aug_logpdf(W[removed_ind, removed_ind], model.nu)
            prior_terms -= normal_logpdf(W[removed_ind, removed_ind], model.w_sigma)

            rpairs = relevant_pairs[removed_ind, removed_ind]
            for p = 1:size(rpairs)[2]
                i = rpairs[1,p]
                j = rpairs[2,p]

                le = latent_effects[i,j]
                oe = observed_effects[i,j]
                old_le = current_latent_effects[i,j]
                for s = 1:length(YY)
                    Y = YY[s]
                    likelihood_terms += log_logit(le + oe, Y[i,j])
                    likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                end
            end 
        elseif L == u
            aug_terms += aug_logpdf(W[new_index,new_index], model.nu)
        else 
            prior_terms += normal_logpdf(W[new_index,new_index], model.w_sigma)
            rpairs = relevant_pairs[new_index, new_index]
            for p = 1:size(rpairs)[2]
                i = rpairs[1,p]
                j = rpairs[2,p]

                le = latent_effects[i,j]
                oe = observed_effects[i,j]
                old_le = current_latent_effects[i,j]
                for s = 1:length(YY)
                    Y = YY[s]
                    likelihood_terms += log_logit(le + oe, Y[i,j])
                    likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                end 
            end 
        end
    else
        for k1 = 1:K # start_index:new_index
            k2_range = k1 != new_index ? new_index :
                       model_spec.symmetric_W ? (1:k1) : 
                       (1:K)
            for k2 = k2_range
                if L <= u
                    aug_terms += aug_logpdf(W[k1,k2], model.nu)
                else
                    prior_terms += normal_logpdf(W[k1,k2], model.w_sigma)
                end
            end
        end

        c_latent_effects = copy(current_latent_effects)

        if L < u || L > u
            assert( L < u || aug_k == new_index)
            for k1 = 1:K
                if model_spec.symmetric_W
                    k2_range = k1 == aug_k ? (1:k1) :
                               k1 < aug_k ? aug_k : []
                else
                    k2_range = k1 == aug_k ? (1:K) : aug_k
                end
                for k2 = k2_range
                    if k1 != new_index && k2 != new_index
                        aug_terms += aug_logpdf(W[k1,k2], model.nu)
                        prior_terms -= normal_logpdf(W[k1,k2], model.w_sigma)
                    end

                    rpairs = relevant_pairs[k1,k2]
                    for p = 1:size(rpairs)[2]
                        i = rpairs[1,p]
                        j = rpairs[2,p]

                        le = latent_effects[i,j]
                        oe = observed_effects[i,j]
                        old_le = c_latent_effects[i,j]
                        # dirty hack to prevent double counting (ie from visiting the same
                        # pair more than once)
                        c_latent_effects[i,j] = le
                        for s = 1:length(YY)
                            Y = YY[s]
                            likelihood_terms += log_logit(le + oe, Y[i,j])
                            likelihood_terms -= log_logit(old_le + oe, Y[i,j])
                        end
                    end
                end
            end
        end
    end

    prior_terms += poisson_logpdf(L,effective_lambda) - 
                   poisson_logpdf(u,effective_lambda)

    (prior_terms, likelihood_terms, aug_terms)
end



function adjust_component_latent_effects(model::ModelState,
                                         model_spec::ModelSpecification,
                                         latent_effects::Array{Float64,2},
                                         relevant_pairs::Array{Array{Int64,2},2},
                                         u::Int,
                                         start_index::Int,
                                         end_index::Int)

    if u < 0
        component_latent_effects = Array(Array{Float64,2}, 1)
        component_latent_effects[1] = latent_effects
        return component_latent_effects
    end

    W = get_W(model, model_spec) 
    (K,K) = size(W)
    (rpK,rpK) = size(relevant_pairs)
#    println("rpK: ", rpK)
#    println("K: ", K)
    assert( rpK == K )

    new_k = end_index + 1
    component_latent_effects = Array(Array{Float64,2}, u+2)


    for u_ind = 1:u+2
        if u_ind <= u
            removed_k = start_index + u_ind - 1
            assert(removed_k <= K)

            component_latent_effects[u_ind] =
                adjust_latent_effects(model, model_spec, latent_effects, relevant_pairs,
                                      u, u-1, removed_k,new_k)
        elseif u_ind == u+1
            component_latent_effects[u_ind] = copy(latent_effects)
        elseif u_ind == u+2
            assert(new_k <= K)
            component_latent_effects[u_ind] =
                adjust_latent_effects(model, model_spec, latent_effects, relevant_pairs,
                                      u, u+1, new_k, new_k)
        end
    end

    component_latent_effects
end

function compute_component_latent_effects(model::ModelState,
                                          model_spec::ModelSpecification,
                                          u::Int,
                                          start_index::Int,
                                          end_index::Int,
                                          node_index::Int)

    if u < 0
        component_latent_effects = Array(Array{Float64,2}, 1)
        component_latent_effects[1] = latent_effects
        return component_latent_effects
    end

    W = get_W(model, model_spec)
    (K,K) = size(W)

    new_k = end_index + 1
    component_latent_effects = Array(Array{Float64,2}, u+2)


    for u_ind = 1:u+2
        if u_ind <= u
            removed_k = start_index + u_ind - 1
            assert(removed_k <= K)

            component_latent_effects[u_ind] =
                compute_latent_effects(model, model_spec, 
                                      u, u-1, removed_k, node_index)
        elseif u_ind == u+1
            component_latent_effects[u_ind] = 
                compute_latent_effects(model, model_spec, 
                                       u, u, 0, node_index)
        elseif u_ind == u+2
            assert(new_k <= K)
            component_latent_effects[u_ind] =
                compute_latent_effects(model, model_spec, 
                                      u, u+1, new_k, node_index)
        end
    end

    component_latent_effects
end


# adjusts latent effects when L != u 
function adjust_latent_effects(model::ModelState,
                               model_spec::ModelSpecification,
                               latent_effects::Array{Float64,2},
                               relevant_pairs::Array{Array{Int64,2},2},
                               u::Int,
                               L::Int,
                               aug_k::Int,
                               new_k::Int)

    W = get_W(model, model_spec)
    (K,K) = size(relevant_pairs)

    assert(aug_k <= K)
    assert(new_k <= K)
    new_latent_effects = copy(latent_effects)
    if L == u - 1
        assert(aug_k < new_k || new_k == 0)

        removed_k = aug_k
        k_range0 = [1:K]
        # latent effects assumes new dimensions haven't been added, so don't remove them!
        if new_k != 0 # only need to worry about new dimension if we have one
            splice!(k_range0, new_k)
        end

        k_range = model_spec.diagonal_W ? removed_k : k_range0

        for k = k_range
            for p = 1:size(relevant_pairs[k, removed_k])[2]
                i = relevant_pairs[k, removed_k][1,p]
                j = relevant_pairs[k, removed_k][2,p]
                new_latent_effects[i,j] -= W[k, removed_k]

            end
            if k != removed_k
                for p = 1:size(relevant_pairs[removed_k, k])[2]
                    i = relevant_pairs[removed_k, k][1,p]
                    j = relevant_pairs[removed_k, k][2,p]
                    new_latent_effects[i,j] -= W[removed_k, k]
                end
            end
        end

    elseif L == u + 1
        new_k = aug_k
        k_range = model_spec.diagonal_W ? new_k : (1:K)
        
        for k = k_range

            if model_spec.symmetric_W
                kmin = min(k,new_k)
                kmax = max(k,new_k)

                if kmax > kmin
                    assert(size(relevant_pairs[kmin,kmax])[2] == 0)
                end
            end

            for p = 1:size(relevant_pairs[k, new_k])[2]
                i = relevant_pairs[k, new_k][1,p]
                j = relevant_pairs[k, new_k][2,p]
                new_latent_effects[i,j] += W[k, new_k]
            end
            if k != new_k
                for p = 1:size(relevant_pairs[new_k, k])[2]
                    i = relevant_pairs[new_k, k][1,p]
                    j = relevant_pairs[new_k, k][2,p]
                    new_latent_effects[i,j] += W[new_k, k]
                end
            end
        end

    end

    new_latent_effects
end

function get_W(model::ModelState,
               model_spec::ModelSpecification)
    if model_spec.positive_W
        return exp(model.weights)
    else
        return copy(model.weights)
    end
end

function get_model_weights(W,
                           model_spec::ModelSpecification)
    if model_spec.positive_W
        return log(W)
    else
        return copy(W)
    end
end

# compute latent effects from scratch
function compute_latent_effects(model::ModelState,
                                model_spec::ModelSpecification,
                                u::Int,
                                L::Int,
                                aug_k::Int,
                                node_index::Int;
                                include_gradient::Bool=false)

    tmodel = copy(model)
    tmodel.tree.nodes[node_index].state = L
    Z = ConstructZ(tmodel.tree)
    w_index_pointers = weight_index_pointers(model.tree)
    start_index = w_index_pointers[node_index]
    new_index = start_index + u
    W = get_W(model, model_spec)
 
    if L < u
        assert( aug_k < new_index)
        W = delete_row_and_col(W, new_index)
        W = delete_row_and_col(W, aug_k)
    elseif L == u
        W = delete_row_and_col(W, new_index)
    elseif L > u
        assert(aug_k == new_index)
    end

    if model_spec.symmetric_W
        W = copy(W)
        symmetrize!(W)
    end

    Z*W*Z'
end

function compute_all_relevant_pairs(model_spec::ModelSpecification,
                                    data::DataState,
                                    K::Int,
                                    Z)

    YY = data.Ytrain
    Ynn = zeros(size(YY[1]))
    Ynn[find(YY[1] .>= 0)] = 1
    Ynn = sparse(Ynn)
    relevant_pairs = [zeros(Int64,(0,0)) for x = 1:K, y = 1:K]

    for k1 = 1:K
        k2_range = model_spec.diagonal_W ? k1 :
                   model_spec.symmetric_W ? (1:k1) :
                   (1:K)
        for k2 = k2_range
            relevant_pairs[k1,k2] = compute_relevant_pairs(model_spec,Ynn,Z,k1,k2)
        end
    end 

    return relevant_pairs
end

# Compute relevant pairs for newly introduced weight parameters
# should be called without changing model.tree from the original
function compute_new_relevant_pairs(model::ModelState,
                                    model_spec::ModelSpecification,
                                    data::DataState,
                                    relevant_pairs::Array{Array{Int64,2},2},
                                    node_index::Int,
                                    start_index::Int,
                                    end_index::Int)
    (K,K) = size(relevant_pairs)
    K += 1
    N::Int = (length(model.tree.nodes)+1)/2
    new_k = end_index + 1

#    println("new_k: ", new_k)
#    println("K: ", K)
    assert(new_k <= K)

    tree = model.tree
    Z = ConstructZ(tree)
    # construct new Z col
    #println("construct Znew")
    leaves = GetLeaves(tree, node_index)
    Zk = zeros(Int64, N)
    Zk[leaves] = 1

    Znew = zeros(Int64, (N, K))
    Znew[:,1:end_index] = Z[:,1:end_index]
    Znew[:,end_index+1] = Zk
    Znew[:,end_index+2:end] = Z[:,end_index+1:end]        
    Znew = sparse(Znew)

    new_relevant_pairs = [zeros(Int64,(0,0)) for x = 1:K, y = 1:K]
    nonzero_element_indices = [x <= end_index ? x : x + 1 for x in 1:K-1]

    YY = data.Ytrain
    Ynn = zeros(size(YY[1]))
    Ynn[find(YY[1] .>= 0)] = 1
    Ynn = sparse(Ynn)

    new_relevant_pairs[nonzero_element_indices, nonzero_element_indices] = copy(relevant_pairs)
    k1_range = model_spec.diagonal_W ? (start_index:end_index+1) : (1:K)
    for k1 = k1_range
        k2 = new_k
        if model_spec.symmetric_W
            k11 = max(k1,k2)
            k22 = min(k1,k2)
            new_relevant_pairs[k11,k22] = compute_relevant_pairs(model_spec,Ynn,Znew,k11,k22)
        else
            new_relevant_pairs[k1,k2] = compute_relevant_pairs(model_spec,Ynn,Znew,k1,k2)
            new_relevant_pairs[k2,k1] = compute_relevant_pairs(model_spec,Ynn,Znew,k2,k1)
        end
    end
    new_relevant_pairs
end

function compute_unaugmented_prob(model::ModelState,
                                  model_spec::ModelSpecification,
                                  relevant_pairs::Array{Array{Int64, 2}, 2},
                                  latent_effects::Array{Array{Float64, 2}, 1},
                                  observed_effects::Array{Float64, 2},
                                  U::Array{Int64, 1},
                                  weight_index_pointers::Array{Int64, 1},
                                  node_index::Int64,
                                  effective_lambda::Float64,
                                  u_index::Int64,
                                  augmented_logprob::Float64)
    W = get_W(model, model_spec) 
    (K,K) = size(W)

    start_index = weight_index_pointers[node_index]
    end_index = start_index + model.tree.nodes[node_index].state - 1
    new_index = end_index + 1

    u = U[u_index]

    current_u = U[end-1]

#    println(U)
#    println(u_index)
    assert(current_u == model.tree.nodes[node_index].state)

    unaugmented_logprob = augmented_logprob
    if model_spec.diagonal_W
        if u < current_u
            removed_ind = start_index + u_index - 1
            unaugmented_logprob -= aug_logpdf(W[removed_ind,removed_ind], model.nu)
            unaugmented_logprob -= aug_logpdf(W[new_index,new_index], model.nu)
        elseif u == current_u 
            unaugmented_logprob -= aug_logpdf(W[new_index,new_index], model.nu)
        end #if u => current_u, no changes are necessary
    else
        if u < current_u
            removed_ind = start_index + u_index - 1
            for k1 = 1:K
                if model_spec.symmetric_W
                    @assert removed_ind < new_index
                    k2_range = k1 == removed_ind || k1 == new_index ? (1:k1) : 
                               k1 < removed_ind ? [] :
                               k1 < new_index ? removed_ind :
                               [removed_ind, new_index]
                    
                else
                    k2_range = k1 == removed_ind || k1 == new_index ? 
                               (1:K) : [removed_ind, new_index]
                end

                for k2 = k2_range
                    unaugmented_logprob -= aug_logpdf(W[k1,k2],model.nu)
                end
            end
        elseif u == current_u
            for k1 = 1:K
                if model_spec.symmetric_W
                    k2_range = k1 == new_index ? (1:k1) : 
                               k1 < new_index ? [] :
                               new_index
                else
                    k2_range = k1 == new_index ? (1:K) : new_index
                end
                for k2 = k2_range
                    unaugmented_logprob -= aug_logpdf(W[k1,k2],model.nu)
                end
            end
        end
    end 

    unaugmented_logprob
end

# Find all pairs i,j such that Z[i,k1]*Z[j,k2] = 1
function compute_relevant_pairs(model_spec::ModelSpecification,
                                mask, #sparse or full binary array, so we only inlude training data
                                Z, #sparse or full binary array
                                k1::Int64,
                                k2::Int64)
    ZZ = Z[:,k1]*Z[:,k2]'
    if model_spec.symmetric_W
        ZZ = ZZ + ZZ' - ZZ .* ZZ'
    end
    #ZZ = ZZ.*mask #only need pairs from Ytrain
    (I, J) = findn(ZZ)
    [I', J']
end

#function symmetrize!(A::Array)
#    triu_inds = find(triu(ones(size(A)),1))
#    A[triu_inds] = tril(A,1)'[triu_inds];
#end

# Doesn't quite handle the symmetric case yet (eg the a[i] and b[j])
function compute_observed_effects(model::ModelState,
                                  model_spec::ModelSpecification,
                                  data::DataState,
                                  i::Int64,
                                  j::Int64)

    X_r = data.X_r
    X_p = data.X_p
    X_c = data.X_c
    observed_effect = model.c

    if model_spec.use_pairwise
        observed_effect += model.beta' * X_r[i,j,:]
    end

    if model_spec.use_parenthood
        if length( X_p[j,:]) > 0
            observed_effect += model.beta_p' * X_p[j,:]
        end
        observed_effect += model.a[j]
    end

    if model_spec.use_childhood
        if length( X_c[i,:]) > 0
            observed_effect += model.beta_c' * X_c[i,:]
        end
        if model_spec.symmetric_W
            observed_effect += model.a[i]
        else
            observed_effect += model.b[i]
        end
    end                                  

    observed_effect
end

function construct_observed_effects(model::ModelState,
                                    model_spec::ModelSpecification,
                                    data::DataState)
    N = size(data.Ytrain[1],1) 
    observed_effects = zeros(Float64, (N,N))
    for i = 1:N
        for j = 1:N
            observed_effects[i,j] = compute_observed_effects(model, model_spec, data, i, j)
        end
    end
    observed_effects
end
# constructs all possible splits of v into 2 sets
function all_splits(v::Array{Int64})
    # Construct power set U, and for each u \in U, construct w = v \ u
    # Then the set of (u,w) will be all possible splits of v

    assert(length(v) < 20)

    U = Vector{Int64}[]
    
    # Construct power set U
    for i = 0:2^(length(v))-1
        x = i
        u = Int64[]
        for j = 1:length(v)
            if bool(x & 1)
                push!(u, v[j])
            end
            x = x >> 1
        end
        push!(U, u)
    end
    (U,reverse(U))
end


#end #profile

