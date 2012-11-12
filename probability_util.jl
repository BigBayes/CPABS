
# Utility probability functions
function normal_logpdf(x)
    -x.*x/2 - log(sqrt(2pi))
end

function poisson_logpdf(k,lambda)
    -lambda + k*log(lambda) - sum(log(1:k))
end

function log_logit(effect, y)
    if y == 1
        return -log(1 + exp(-effect))
    elseif y == 0
        return -effect - log(1 + exp(-effect))
    end
    return 0.0
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
