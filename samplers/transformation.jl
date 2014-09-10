
# Transformation type for handling transformations of variables
# Note that the parent MCMCState is included as an argument for
# transformed_gradient and transformation_logjacobian only so that
# the Jacobian for the transformed parameter may depend on other parameters  
type Transformation
    # transformation mapping from represented space to original space
    transformation_function::Function
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

function identity_transformed_gradient(th,g)
    g
end

function identity_transformation_logjacobian(th)
    0.0
end

IdentityTransformation = Transformation(identity_transformation,
                                        identity_transformed_gradient,
                                        identity_transformation_logjacobian)


function logspace_transformation(th)
    exp(th)
end

function logspace_transformed_gradient(h,g)
    exp(th).*g + ones(size(th))
end

function logspace_transformation_logjacobian(th)
    sum(th)
end

LogSpaceTransformation = Transformation(logspace_transformation,
                                        logspace_transformed_gradient,
                                        logspace_transformation_logjacobian)


# Log Transformation for only diagonal entries in a matrix
function diag_logspace_transformation(th)
    N2 = length(th)
    N = int(sqrt(N2))
    ret = deepcopy(th)
    dind = diagind(N,N)
    ret[dind] = exp(th[dind])
    ret
end

function diag_logspace_transformed_gradient(th,g)
    N2 = length(th)
    N = int(sqrt(N2))
    ret = deepcopy(g)
    dind = diagind(N,N)
    ret[dind] = exp(th[dind]).*g[dind] + ones(length(th[dind]))
    ret
end

function diag_logspace_transformation_logjacobian(th)
    N2 = length(th)
    N = int(sqrt(N2))
    dind = diagind(N,N)
    sum(th[dind])
end

DiagonalLogSpaceTransformation = Transformation(diag_logspace_transformation,
                                                diag_logspace_transformed_gradient,
                                                diag_logspace_transformation_logjacobian)

# Reduced natural transformation for parameters with domain [0,1]
function reduced_natural_transformation(th)
    # different calculations depending on sign for numerical stability
    if th > 0.0
        1.0./(1.0.+exp(-th))
    else
        exp(th)./(1.0.+exp(th))
    end
end

function reduced_natural_transformed_gradient(th,g)
    if th > 0.0
        exp(-th)./(1.+exp(-th))^2.*g .+ 1.0 .- 2.0./(1.+exp(-th))
    else
        exp(th)./(1.+exp(th))^2.*g .+ 1.0 .- 2.0exp(th)./(1.+exp(th))
    end
end

function reduced_natrual_transformation_logjacobian(th)
    if th > 0.0
        -th .+ 2.0log(1.+exp(-th))
    else
        th .+ 2.0log(1.+exp(th))
    end
end

ReducedNaturalTransformation = Transformation(reduced_natural_transformation,
                                              reduced_natural_transformed_gradient,
                                              reduced_natrual_transformation_logjacobian)
