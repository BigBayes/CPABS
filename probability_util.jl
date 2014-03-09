require("model.jl")

function incomplete_gamma(t::Float64, b::Float64)
    quadgk(x -> x^(t-1)*exp(-x), zero(t), b, reltol=100*eps(t))[1]
end


# Utility probability functions
function normal_logpdf(x, sigma)
    -x.*x/(2sigma*sigma) - 0.5log(2pi) - log(sigma)
end
function normal_logpdf_dx(x, sigma)
    -x/(sigma*sigma)
end

function normal_rand(scale)
    rand(Normal(0.0, scale))
end

function exp_logpdf(x, scale)
    -x/scale - log(scale)
end

function exp_logpdf_dx(x, scale)
    -ones(size(x))/scale
end

function exp_rand(scale)
    rand(Exponential(scale))
end

function poisson_logpdf(k,lambda)
    if lambda == 0.0 && k == 0
        return 0.0
    end
    -lambda + k*log(lambda) - sum(log(1:k))
end

# Must be normalized for use in expanding dimensions for variable dimension slice sampling
function aug_logpdf(x, nu)
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

function log_logistic_dx(effect, y)
    -(1-y) + exp(-effect)/(1+exp(-effect))
end


function log_predictive(effect)
    -log(1+exp(-effect))
end

# Probability Manipulation

function logsumexp(x)
    max_x = max(x)
    xp = x - max_x
    log(sum(exp(xp))) + max_x
end

function logsumexp_d_dx(x,xp)
    #zeros in shape of xp[1]
    xout = xp[1]-xp[1]

    max_x = max(x)
    xs = x - max_x

    exp_xs = exp(xs)
    sum_exp_xs = sum(exp_xs)
    for i = 1:length(xp)
        xout += exp_xs[i]*xp[i]
    end

    xout /= sum_exp_xs

    #sum(exp(xs).*xp) / sum(exp(xs))
    xout
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

# Prediction

# test metrics.  AUC computation trick found in Konstantina Palla's ILA code
function error_and_auc(logit_args::Array{ Array{Float64, 1}, 2},
                       data::DataState)

    YY = data.Ytrain
    (N,N) = size(YY[1])

#    minargs = min( [ min(logit_args[i,j]) for i = 1:N, j = 1:N])
#    maxargs = max( [ max(logit_args[i,j]) for i = 1:N, j = 1:N])
#
#    bias_range = linspace(-maxargs-2,-minargs+2,100)
#
#
#    roc_false_positives = zeros(Int, length(bias_range))
#    roc_true_positives = zeros(Int, length(bias_range))
#    roc_positives = zeros(Int, length(bias_range))
#
    probs = zeros(N,N)
    train01error = 0
    test01error = 0

    for i = 1:N
        for j = 1:N
            logit_arg_list = logit_args[i,j]
            (prob, logprob) = averaged_prediction(logit_arg_list)
            probs[i,j] = prob

            for s = 1:length(YY)
                Ytrain = YY[s]
                Ytest = data.Ytest[s]
                if Ytrain[i,j] >= 0
                    train01error += round(prob) != Ytrain[i,j]
                elseif Ytest[i,j] >= 0
                    test01error += round(prob) != Ytest[i,j]
                end

            end
        end
    end

    Ytrain_flat = Int64[]
    Ytest_flat = Int64[]
    probs_flat = Float64[]

    for s = 1:length(YY)
        append!(Ytrain_flat, YY[s][:])
        append!(Ytest_flat, data.Ytest[s][:])
        append!(probs_flat, probs[:])
    end

    Itrain = find(Ytrain_flat .>= 0)
    Itest = find(Ytest_flat .>= 0)
    test_pos = find(Ytest_flat .== 1)
    test_neg = find(Ytest_flat .== 0)

    test_probs = probs_flat[Itest]
    p = sortperm(test_probs)
    Ytest_sorted = Ytest_flat[Itest[p]]

    num_links = sum(Ytest_sorted .== 1)
    num_nonlinks = sum(Ytest_sorted .== 0)

    # sum of ranks (minus correction) gives the true positive count, false positive count area
    count_under_curve = sum(find(Ytest_sorted .== 1)) - (num_links*(num_links+1))/2

    auc = (count_under_curve ) / (num_links * num_nonlinks)

    train_error_rate = train01error / length(Itrain)
    test_error_rate = test01error / length(Itest)


    (train_error_rate, test_error_rate, auc)
end

function averaged_prediction(logit_arg::Array{Float64, 1})

    logprobs = log_predictive(logit_arg)
    logprob = logsumexp(logprobs) - log(length(logprobs))
    (exp(logprob), logprob)
end



function trapz(x::Array{Float64,1}, y::Array{Float64,1})
    assert(length(x) == length(y))

    inds = find(!isnan(x) .* !isinf(x) .* !isnan(y) .* !isinf(y))

    xx = x[inds]
    yy = y[inds]

    ly = yy[1:end-1]
    ry = yy[2:end]
   
    lx = xx[1:end-1]
    rx = xx[2:end] 

    sum(0.5*(ry + ly) .* (rx-lx))
end

