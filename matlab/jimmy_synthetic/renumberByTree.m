function [newNumbering Znew, Xnew] = renumberByTree(Z, X)
    %assumes that the features are generated in order of recursing the
    %tree, i.e. if the previous feature is not there then it is not a
    %descendant of a feature already processed
    
    Znew = Z;
    numFeatures = size(Z,2);
    newNumbering = 1:size(Z,1);
    lastFeatureTop = 1;
    afterLastFeature = -1; %uninitialized
    for k = 1:numFeatures;
        onFeatures = find(Znew(:,k));
        featureStart = -1;
        if (k == 1 || Znew(onFeatures(1), k-1) == 1)
            %put it at the top of the current Feature
            featureStart = lastFeatureTop;
        else
            %put it after the current feature
            featureStart = afterLastFeature;
        end
        ZtoMove = Znew(onFeatures, :);
        Znew(onFeatures,:) = [];
        Znew = [Znew(1:featureStart-1,:); ZtoMove; Znew(featureStart:end,:)];
        numbersToMove = newNumbering(onFeatures);
        newNumbering(onFeatures) = [];
        newNumbering = [newNumbering(1:featureStart-1) numbersToMove newNumbering(featureStart:end)];
        lastFeatureTop = featureStart;
        afterLastFeature = featureStart + length(onFeatures);
    end
    Xnew = zeros(size(X));
    numActors = size(X,1);
    for i = 1:numActors
        for j = 1:numActors
            %Xnew(newNumbering(i), newNumbering(j)) = X(i,j);
            Xnew(i,j) = X(newNumbering(i), newNumbering(j));
        end
    end
    assert(nnz(Xnew) == nnz(X));
    W = normrnd(0,1, numFeatures, numFeatures);
    ll_old = sum(sum(X .*log(1./(1 + exp(-Z * W * Z'))) + (1-X).*log(1-1./(1 + exp(-Z * W * Z'))) ));
    ll_new = sum(sum(Xnew .*log(1./(1 + exp(-Znew * W * Znew'))) + (1-Xnew).*log(1-1./(1 + exp(-Znew * W * Znew'))) ));
    assert(abs(ll_old - ll_new) < 0.00001);
end


% function newNumbering = renumberByTree(children, currentNode, newNumbering)
%     if children{currentNode}(1) == -1; %leaf node
%         %node indices for leaves correspond to actor indices, since they are placed into the tree first
%         newNumbering = [newNumbering currentNode];
%         return;
%     end
% 
%     for a = 1:length(children{currentNode})
%         newNumbering = renumberByTree(children, children{currentNode}(a), newNumbering);
%     end
% end