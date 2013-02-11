%% function for reading in .jla files.
% metrics files:  These contain train/test errors, log likelihoods, AUCs, numbers of 
% feautes (ie Ks).  Specifically from one run we have in this order:
% iterations:  These are the iterations on which the model was recorded (trees,features,
%              and weights), and the more expensive metrics were recored if past burn in
%              (ie AUCs, errors, etc).
% train_errors: averaged over all (non-burnin) available prior samples
% test_errors
% averaged_test_LL: predicted test log likelihood averaged over all samples up to now
% AUC: roc AUC over time
% K: number of features
% train_LL: train log likelihood for this iteration.  available for all iterations
% test_LL: test log likelihood for this iteration.

% trees files:  The trees from the run, one per element in metrics->iterations
% Format: N-1 x 4 matrix (same as is used in Python)
% first two columns are indices of children joined to make the current node
% third column is the time
% fourth column is the total number of descendants

% weights files: The weight matrices, one per element in metrics->iterations
% features files: The feature counts per node of the tree.
% data files: The train test splits, for this run.
%             (-1s correspond to missing, 0/1 correspond to no-link/link)
function [arrays] = restore(filename)
fid = fopen(filename,'r');
array_size = fread(fid, 1, 'int64');
nds = fread(fid, array_size, 'int64');
lengths = fread(fid, sum(nds), 'int64');

sizes = cell(array_size,1);

count = 0;
for i = 1:array_size
    sz = zeros( 0,0, 'int64');
    for j = 1:nds(i)
        count = count + 1;
        sz = [sz; lengths(count)];
    end
    sizes{i} = sz;
end

arrays = cell(array_size,1);
for i = 1:array_size
    arrays{i} = fread(fid, sizes{i}', 'double');
end

fclose(fid);

