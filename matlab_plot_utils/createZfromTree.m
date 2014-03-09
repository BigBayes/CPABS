function Z = createZfromTree(tree, features)
    %go backwards from the tree to a Z matrix.  The tree is in Levi's
    %format, not Jimmy's one used to create synthetic trees (which are
    %roughly equivalent formats).
    %@author Jimmy, with debugging help from Levi
    numActors = size(tree,1) + 1;
    Z = zeros(numActors, sum(features));
    Z = addFeaturesBelow(tree, features, Z, 0, [], 2 .* numActors - 1);
end

function [Z maxFeature] = addFeaturesBelow(tree, features, Z, maxFeature, featuresAboveCurrentNode, currentNode)
%recursively pass through the tree and create features
    numActors = size(Z,1);

    %add features on this branch
    for i = 1:features(currentNode)
        featuresAboveCurrentNode = [featuresAboveCurrentNode maxFeature + 1];
        maxFeature = maxFeature + 1;
        if maxFeature == 2
            1+1;
        end
    end
    
    if currentNode <= numActors %leaf node
        for i = 1:length(featuresAboveCurrentNode)
            Z(currentNode, featuresAboveCurrentNode(i)) = 1;
        end
        return;
    end
    [Z maxFeature] = addFeaturesBelow(tree, features, Z, maxFeature, featuresAboveCurrentNode, round(tree(currentNode-numActors, 1))); %left
    [Z maxFeature] = addFeaturesBelow(tree, features, Z, maxFeature, featuresAboveCurrentNode, round(tree(currentNode-numActors, 2))); %right
    
    
end