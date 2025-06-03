%% conduct a single round simulation for part of recordings
function [simulation_result, real_result] = single_round_simulation_part_recordings(...
    theta, other_settings, data_path_collections, fps_org, fps_upsampled, ~, ~)
% data_path_collections is either trainSet or testSet from
% createStratifiedTrainTestSplits, indexed by cv_splits(i).train or
% cv_splits(i).test
numFiles = numel(data_path_collections);

% — preallocate sliced outputs —
simulation_result_tmp = cell(1, numFiles);
real_result_tmp       = cell(1, numFiles);

% — parallel loop —
parfor ii = 1:numFiles
    % 1) load + process behavioral data
    aqua_data_path = data_path_collections{ii};
    [resp_seq, run_seq, stim_onset, reward_onset, TT] = ...
        behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled);

    % 2) precompute inputs
    [pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun] = ...
        precompute_inputs(theta, other_settings, aqua_data_path, fps_org, fps_upsampled);

    % 3) run simulation
    [chemical_table, ~] = single_round_simulation_one_recording(...
        theta, other_settings, run_seq, stim_onset, reward_onset, TT, ...
        pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun);
    c_cyto = chemical_table{:,"cytosolic_Ca"};

    % 4) store the sliced outputs
    simulation_result_tmp{ii} = c_cyto;
    real_result_tmp{ii}       = resp_seq;
end

% — post‑process: turn cell arrays into structs with nice fieldnames —
simulation_result = struct();
real_result       = struct();
for ii = 1:numFiles
    path = data_path_collections{ii};
    [nm, loc, num] = extractInfoFromPath(path);
    field = sprintf('%s_%s_%s', nm, loc, num);
    simulation_result.(field) = simulation_result_tmp{ii};
    real_result.(field)       = real_result_tmp{ii};
end

end