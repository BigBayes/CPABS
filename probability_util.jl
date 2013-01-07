
# Utility probability functions
function normal_logpdf(x, sigma)
    -x.*x/(2sigma*sigma) - 0.5log(2pi) - log(sigma)
end

function poisson_logpdf(k,lambda)
    -lambda + k*log(lambda) - sum(log(1:k))
end

# Must be normalized for use in expanding dimensions for variable dimension slice sampling
function t_logpdf(x, nu)
    -x.*x/(2nu*nu) - 0.5log(2pi) - log(nu)
    #lgamma(0.5*(nu+1)) - lgamma(0.5nu) - 0.5*log(nu*pi) - 0.5*(nu+1)*log(1+x*x/nu)
end

# Sigma = I
function multivariate_t_logpdf(x_squared_norm, p, nu)
   lgamma(0.5*(nu+p)) - lgamma(0.5nu) - 0.5p*log(nu*pi) - 0.5*(nu+p)*log(1 + x_squared_norm/nu)
end

function log_logit(effect, y)
    value = 0.0
    if y == 1
        value = -log(1 + exp(-effect))
        if isinf(value)
            value = effect
        end
    elseif y == 0
        value = -effect - log(1 + exp(-effect))
        if isinf(value)
            value = 0.0
        end
    end
    return value
end


# Probability Manipulation

function logsumexp(x)
    max_x = max(x)
    xp = x - max_x
    log(sum(exp(xp))) + max_x
end

function exp_normalize(x)
    xp = x - max(x)
    exp_x = exp(xp)
    exp_x / sum(exp_x)
end

# Basic Sampling

function randmult(x)
    v = cumsum(x)
    assert( abs(v[end] - 1.0) < 10.0^-8)
#    if v[end] != 1.0
#        println("v[end]: ", v[end])
#        assert(v[end] == 1.0)
#    end

    u = rand()
    i = 1
    while u > v[i]
        i += 1
    end
    i
end

function randpois(lambda)
    L = exp(-lambda)
    k = 0
    p = 1.0
    while p > L
        k += 1
        p = p * rand()
    end
    k - 1
end
