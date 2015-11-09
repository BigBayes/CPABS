
# Transformation type for handling transformations of variables
# Note that the parent MCMCState is included as an argument for
# transformed_gradient and transformation_logjacobian only so that
# the Jacobian for the transformed parameter may depend on other parameters  
type Transformation
    # transformation mapping from represented space to original space
    transformation_function::Function
    # inverse transformation mapping from original space to represented space
    transformation_inverse::Function
    # Maps gradient of logdensity computed in original space to transformed space,
    # including the log Jacobian term.  This function takes the transformed (ie represented)
    # variable as input, not the variable in the original space. 
    transformed_gradient::Function
    # log of the Jacobian determinant of the transformation
    transformation_logjacobian::Function
end


function identity_transformation(th)
    th
end

function identity_inverse(x)
    x
end

function identity_transformed_gradient(th,g)
    g
end

function identity_transformation_logjacobian(th)
    0.0
end

IdentityTransformation = Transformation(identity_transformation,
                                        identity_inverse,
                                        identity_transformed_gradient,
                                        identity_transformation_logjacobian)


function logspace_transformation(th)
    exp(th)
end

function logspace_inverse(x)
    log(x)
end

function logspace_transformed_gradient(th,g)
    exp(th).*g + ones(size(th))
end

function logspace_transformation_logjacobian(th)
    abs(sum(th))
end

LogSpaceTransformation = Transformation(logspace_transformation,
                                        logspace_inverse,
                                        logspace_transformed_gradient,
                                        logspace_transformation_logjacobian)



# Reduced natural transformation for parameters with domain [0,1]
function reduced_natural_transformation(th)
    # different calculations depending on sign for numerical stability
    result = zeros(size(th))
    for t = 1:length(th)
        if th[t] > 0.0
            result[t] = 1.0/(1.0+exp(-th[t]))
        else
            result[t] = exp(th[t])/(1.0+exp(th[t]))
        end
    end
    result
end

function reduced_natural_inverse(x)
    result = zeros(size(x))
    for t = 1:length(x)
        if x[t] > 0.5
            result[t] = -log(1.0/x[t] - 1.0)
        else
            result[t] = log(x[t]/(1.0-x[t]))
        end
    end
    result
end

function reduced_natural_transformed_gradient(th,g)
    result = zeros(size(th))
    logjacobian = reduced_natrual_transformation_logjacobian(th, take_absolute_value=false)
    sgn = sign(logjacobian)
    for t = 1:length(th)
        if th[t] > 0.0
            result[t] = exp(-th[t])/(1+exp(-th[t]))^2*g[t] + sgn*(1.0 - 2.0/(1+exp(-th[t])))
        else
            result[t] = exp(th[t])/(1+exp(th[t]))^2*g[t] + sgn*(1.0 - 2.0exp(th[t])/(1+exp(th[t])))
        end
    end
    result
end

function reduced_natrual_transformation_logjacobian(th; take_absolute_value = true)
    result = 0.0
    for t = 1:length(th)
        if th[t] > 0.0
            result += -th[t] + 2.0log(1+exp(-th[t]))
        else
            result += th[t] + 2.0log(1+exp(th[t]))
        end
    end
    if take_absolute_value
        return abs(result)
    else
        return result
    end
end

ReducedNaturalTransformation = Transformation(reduced_natural_transformation,
                                              reduced_natural_inverse,
                                              reduced_natural_transformed_gradient,
                                              reduced_natrual_transformation_logjacobian)
