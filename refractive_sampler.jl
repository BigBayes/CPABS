require("Options")
using OptionsMod

function refractive_sampler(init_x::Vector{Float64},
                            f::Function,
                            grad_f::Function,
                            opts::Options) 
  @defaults opts w=.01 refractive_index_ratio=nothing m=10 sample_m=false w_dist=nothing 
  

  if w_dist != nothing
    w = w* w_dist()
  end

#  if refractive_index_ratio < 1.0
#    error("refractive_index_ratio should be >= 1.0")
#  end

  if w <= 0
    error("Negative w not allowed")
  end
  if m <= 0
    error("Limit on steps must be positive")
  end

  if sample_m
    m = rand(Poisson(m))
  end


  local rir = refractive_index_ratio
  local accept_logprob = 0.0
  local d
  local x0, x
  local p

  local cos_th_1, cos_th_2, cos2_th_2
  local u, r1, r2
  local grad, grad_fx
  local fx, fx0


  local jacobian_term = 0.0


  x0 = copy(init_x)
  x = copy(x0)

  d = length(x0)
  p = randn(d)


  fx0= f(x0)
  accept_logprob = 0.0

  for iteration = 1:m+1
    fx = f(x)
    if !isfinite(fx)
        println("reject: out of bounds?")
        accept_logprob = -Inf
        break
    end
    grad = grad_f(x)


    grad_norm = sqrt(dot(grad,grad))
    grad = grad/grad_norm

    Cg = grad_norm*w
    if dot(p, grad) > 0.0
        u = grad
        cos_th_1 = dot(p,u)/sqrt(dot(p,p))
        Ct2 = cos_th_1^2
       
        r = rir == nothing ? 1 + 4*Cg*Ct2/(4*Ct2*(d-1) + 1-Ct2) : rir
        r1 = 1.0
        r2 = r #rir(grad_norm, w, d)
    else
        u = -grad
        cos_th_1 = dot(p,u)/sqrt(dot(p,p)) 
        Ct2 = cos_th_1^2

        r = rir == nothing ? 1 + 4*Cg*Ct2/(4*Ct2*(d-1) + 1-Ct2) : rir
        r1 = r #rir(grad_norm, w, d)
        r2 = 1.0
    end

    cos2_th_2 = 1 - (r1/r2)^2*(1.0-cos_th_1^2)

    jacobian_term = 0.0
    if cos2_th_2 < 0.0 # reflect
        p = p - 2dot(p,u) * u
    else # refract
        cos_th_2 = sqrt(cos2_th_2)
        p = r1/r2*p - sqrt(dot(p,p))*(r1/r2* cos_th_1 - cos_th_2)*u
        jacobian_term = (d-1)*(log(r1) - log(r2)) + log(cos_th_1) - log(cos_th_2)
    end

    accept_logprob += jacobian_term
    if iteration == m+1
        accept_logprob += fx - fx0
    else
        x += w*p
    end

  end


#  println("fx: $fx")
#  println("fx0: $fx0")
#  println("accept_logprob: $accept_logprob")

  if rand() < exp(accept_logprob)
    return x
  else
    return x0
  end


end
