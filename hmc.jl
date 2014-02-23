require("Options")
using OptionsMod

function hmc_sampler(x::Vector{Float64},
                     U::Function,
                     grad_U::Function,
                     opts::Options)


  @defaults opts stepsize=.01 numsteps=50 nadapt=0 w=.01 

  current_U = U(x)

  local epsilon = stepsize
  local L = numsteps


  current_q = copy(x)
  q = copy(current_q)

  p = randn(length(q))  # independent standard normal variates
  current_p = deepcopy(p)

  gradient = grad_U(x)
  # Make a half step for momentum at the beginning
  p = p + epsilon * gradient / 2

  # Alternate full steps for position and momentum
  for i in 1:L
    # Make a full step for the position

    q = q + epsilon * p

    Ux = U(q)
    if !isfinite(Ux)
        break
    end

    # Make a full step for the momentum, except at end of trajectory
    if i!=L
        gradient = grad_U(q) 
        p = p + epsilon * gradient
    end

  end

  # Make a half step for momentum at the end.
  proposed_U = U(q)
  gradient = grad_U(q) 

  p = p + epsilon * gradient / 2

  # Negate momentum at end of trajectory to make the proposal symmetric
  # p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory

  current_K = 0.5 * dot(current_p,current_p) 
  proposed_K = 0.5 * dot(p,p) 
  H_initial = -current_U + current_K
  H_proposed = -proposed_U + proposed_K

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept_prob = exp(-(H_proposed - H_initial))
  if isnan(accept_prob)
    println("HMC")
    println("H_prop: $H_proposed")
    println("H_init: $H_initial")
    println("current U, K: $current_U, $current_K")
    println("proposed U, K: $proposed_U, $proposed_K")
    println("gradient_norm: $(sqrt(dot(gradient,gradient) ))")
  end


  if rand() < accept_prob 
    return q 
  else
    return current_q
  end



end

