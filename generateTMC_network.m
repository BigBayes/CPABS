function [X, W, Z, parents, children, featuresOnBranchAboveNode] = generateTMC_network(branchLengthGamma, numActors, featureRateLambda, LFRM_sigma, interceptEpsilon)
%generate from the time marginalized coalescent for networks model.
%@author Jimmy Foulds
%X: the network generated
%W: LFRM weight matrix, numFeatures x numFeatures
%Z: numActors x numFeatures latent feature matrix
%parents: an array, for every node in the TMC tree, its index of its parent
%children: a cell array, for every node in the TMC tree i, children{1} is its left child and children{2} is its right child.
%featuresOnBranchAboveNode: A cell array, for every node in the TMC tree, the indices of features that appeared in the branch above it.

    [parents, children] = generateTMC_tree(numActors);
    [Z, featuresOnBranchAboveNode] = generateTMC_features(numActors, parents, children, branchLengthGamma, featureRateLambda);
    numFeatures = size(Z,2);
    W = normrnd(0,LFRM_sigma,numFeatures,numFeatures);
    X = rand(numActors) < 1./(1 + exp(-(Z * W * Z' + interceptEpsilon)));
end

function [parents, children] = generateTMC_tree(numActors)
%Make the tree itself via the TMC (i.e. joining nodes randomly until we are done)
    children = cell(numActors,1);
    for i = 1:numActors
        children{i} = [-1,-1];
    end
    openList = 1:numActors;
    while length(openList) > 1
        l_ind = randi(length(openList));
        l = openList(l_ind);
        r_ind = l_ind;
        while r_ind == l_ind
            r_ind = randi(length(openList));
        end
        r = openList(r_ind);
        newNode = length(children) + 1;
        children{newNode} = [l r];
        parents(newNode) = -1;
        parents(l) = newNode;
        parents(r) = newNode;
        
        openList([l_ind r_ind]) = [];
        openList(end + 1) = newNode;
    end
end

function [Z, featuresOnBranchAboveNode] = generateTMC_features(numActors, parents, children, branchLengthGamma, featureRateLambda)
%Generate the feature matrix.  Nice wrapper function for the recursive function.
    [Zsparse_i, Zsparse_j, featuresOnBranchAboveNode, maxFeature] = generateFeaturesBelow(numActors, parents, children, branchLengthGamma, featureRateLambda, cell(length(children),1), [], [], 0, [], length(children), 0, 1);
    Z = sparse(Zsparse_i, Zsparse_j, ones(length(Zsparse_i),1), numActors, maxFeature);
end


function [Zsparse_i, Zsparse_j, featuresOnBranchAboveNode, maxFeature] = generateFeaturesBelow(numActors, parents, children, branchLengthGamma, featureRateLambda, featuresOnBranchAboveNode, Zsparse_i, Zsparse_j, maxFeature, featuresFromAbove, currentNode, currentTime, numAncestors)
%recursive function for the tree, generate the features in the portion of the tree starting from currentNode
    if children{currentNode}(1) == -1; %leaf node
        branchLength = 1 - currentTime;
        numFeaturesAboveNode = poissrnd(featureRateLambda .* branchLength);
        featuresOnBranchAboveNode{currentNode} = (maxFeature+1:maxFeature+numFeaturesAboveNode)';
        
        Zsparse_i = [Zsparse_i; ones(numFeaturesAboveNode + length(featuresFromAbove), 1) .* currentNode]; %node indices for leaves correspond to actor indices, since they are placed into the tree first
        Zsparse_j = [Zsparse_j; featuresFromAbove; (maxFeature+1:maxFeature+numFeaturesAboveNode)'];
        maxFeature = maxFeature + numFeaturesAboveNode;
        return;
    end
    branchLength = (1-branchLengthGamma).* branchLengthGamma^(numAncestors-1);
    numFeaturesAboveNode = poissrnd(featureRateLambda .* branchLength);
    featuresOnBranchAboveNode{currentNode} = (maxFeature+1:maxFeature+numFeaturesAboveNode)';

    for a = 1:length(children{currentNode})
        [Zsparse_i, Zsparse_j, featuresOnBranchAboveNode, maxFeature] = generateFeaturesBelow(numActors, parents, children, branchLengthGamma, featureRateLambda, featuresOnBranchAboveNode, Zsparse_i, Zsparse_j, maxFeature + numFeaturesAboveNode, [featuresFromAbove; featuresOnBranchAboveNode{currentNode}], children{currentNode}(a), currentTime + branchLength, numAncestors + 1);
    end
end