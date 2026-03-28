%% Evaluate a population of candidates using flattened parallelism
% Instead of parfeval per candidate (each with inner parfor over recordings),
% this flattens to parfeval per (candidate, recording) pair.
% All pop_size * numFiles tasks run at the same level, maximizing core usage.
%
% Inputs:
%   candidates       - cell array of containers.Map (1 x pop_size)
%   other_settings   - settings map
%   data_paths       - cell array of recording paths (numFiles x 1)
%   sliced_behavior  - cell array of behavior structs (numFiles x 1)
%   fps_org          - original fps
%   fps_upsampled    - upsampled fps
%   save_path        - path for saving intermediate results
%
% Output:
%   cachedFitness    - cell array (1 x pop_size) of metric column vectors
function cachedFitness = evaluate_population_flat(candidates, other_settings, ...
    data_paths, sliced_behavior, field_names, fps_org, fps_upsampled, save_path)

pop_size = numel(candidates);
numFiles = numel(data_paths);
totalTasks = pop_size * numFiles;

pool = gcp();

% Submit all (candidate, recording) pairs as individual tasks
futures = parallel.FevalFuture.empty(totalTasks, 0);
for cc = 1:pop_size
    for rr = 1:numFiles
        idx = (cc - 1) * numFiles + rr;
        futures(idx) = parfeval(pool, @simulate_one_task, 1, ...
            candidates{cc}, other_settings, data_paths{rr}, ...
            sliced_behavior{rr}, fps_org, fps_upsampled);
    end
end

% Gather all results
sim_results = cell(totalTasks, 1);
for idx = 1:totalTasks
    sim_results{idx} = fetchOutputs(futures(idx));
end

% Aggregate per candidate: build sim_result and real_result structs,
% then compute multi-objective metrics
cachedFitness = cell(1, pop_size);
for cc = 1:pop_size
    sim_struct  = struct();
    real_struct = struct();
    for rr = 1:numFiles
        idx = (cc - 1) * numFiles + rr;
        fn = field_names{rr};
        sim_struct.(fn)  = sim_results{idx};
        real_struct.(fn) = sliced_behavior{rr}.resp_seq;
    end
    cachedFitness{cc} = multi_metric_eval(sim_struct, real_struct, save_path);
end

end
