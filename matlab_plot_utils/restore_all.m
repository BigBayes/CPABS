function [results] = restore_all(result_path, id_string, num_trials)

    results = struct;
    results.metrics = metrics_array2struct(...
        restore_in_range([result_path,'/metrics_',id_string], 1:num_trials)); 
    results.trees = restore_in_range([result_path,'/trees_',id_string], 1:num_trials); 
    results.features = restore_in_range([result_path,'/features_',id_string],1:num_trials); 
    results.weights = restore_in_range([result_path,'/weights_',id_string], 1:num_trials); 
    results.datas = restore_in_range([result_path,'/data_',id_string], 1:num_trials); 

end

function [results] = restore_in_range(filename_base, index_range)
    results = cell(length(index_range),1);
    for i = index_range
        results{i} = restore([filename_base,num2str(i),'.jla']);
    end
end

function [s]  = metrics_array2struct(metrics)
    s = cell(length(metrics),1);

    for i = 1:length(metrics)
        s{i} = struct;
        metric = metrics{i};
        s{i}.iters = metric{1};
        s{i}.train_errors = metric{2};
        s{i}.test_errors = metric{3};
        s{i}.avg_test_LLs = metric{4};
        s{i}.AUCs = metric{5};
        s{i}.Ks = metric{6};
        s{i}.trainLLs = metric{7};
        s{i}.testLLs = metric{8};
    end
end
