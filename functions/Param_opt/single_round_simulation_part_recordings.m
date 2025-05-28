%% conduct a single round simulation for part of recordings
function [simulation_result, real_result] = single_round_simulation_part_recordings(...
    theta, other_settings, data_path_collections, fps_org, fps_upsampled, save_tables, important_fields)
% data_path_collections is either trainSet or testSet from
% createStratifiedTrainTestSplits, indexed by cv_splits(i).train or
% cv_splits(i).test
numFiles = numel(data_path_collections);

% — preallocate sliced outputs —
simulation_result_tmp = cell(1, numFiles);
real_result_tmp       = cell(1, numFiles);

% — base folder for saving tables —
arch_type = computer('arch');
if strcmp(arch_type, 'win64')
    baseSaveDir = 'C:\Users\james\CBIL\Astrocyte\Temporal_model_eval';
elseif strcmp(arch_type, 'glnxa64')
    baseSaveDir = '/work/Jamespeng/Astrocyte/Temporal_model_eval';
end

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
    [chemical_table, flux_table] = single_round_simulation_one_recording(...
        theta, other_settings, run_seq, stim_onset, reward_onset, TT, ...
        pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun);
    c_cyto = chemical_table{:,"cytosolic_Ca"};

    % 4) store the sliced outputs
    simulation_result_tmp{ii} = c_cyto;
    real_result_tmp{ii}       = resp_seq;

    % 5) compute scalar correlation
    corr_val = corr(c_cyto(:), resp_seq(:));  % single r

    if save_tables
        S = struct()
        % — parse out mouse, RecLoc, recording name —
        [parentFolder, fileBase, ~] = fileparts(aqua_data_path);
        parts      = strsplit(parentFolder, filesep);
        mouse_name = parts{end-1};
        RecLoc     = parts{end};
        tb         = strsplit(fileBase, '_');
        recording  = tb{1};

        % — ensure output folder exists —
        outDir = fullfile(baseSaveDir, mouse_name, RecLoc, recording);
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        % — pack into one struct —
        S.chemical_table = chemical_table;
        S.flux_table     = flux_table;
        S.resp_seq       = resp_seq;
        S.run_seq        = run_seq;
        S.stim_onset     = stim_onset;
        S.reward_onset   = reward_onset;
        S.TT             = TT;

        % — build filename —
        if important_fields
            fname = sprintf('GA_train_test_important_fields_%.3f.mat', corr_val);
        else
            fname = sprintf('GA_train_test_all_fields_%.3f.mat', corr_val);
        end
        fullFname = fullfile(outDir, fname);

        % — call our helper (legal inside parfor) —
        parsaveStruct(fullFname, S);
    end
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

%% helper function
function [extractedName, extractedLocation, extractedNumber] = extractInfoFromPath(path)
% extractInfoFromPath Extracts the name, location, and number from a given file path.

% Split the file path into parts using the system-specific file separator.
parts = strsplit(path, filesep);
% Locate the directory 'AQUA_processed_data'
idx = find(strcmp(parts, 'AQUA_processed_data'), 1);
% Check that the expected parts exist in the path
if isempty(idx) || numel(parts) < idx + 3
    error('The given path does not follow the expected format.');
end
% Extract the name (folder immediately after 'AQUA_processed_data')
extractedName = parts{idx+1};
% Extract the location (the folder after the name)
extractedLocation = parts{idx+2};
% Extract the file name (should be the next element in the parts)
% Remove the extension using fileparts
[~, fileBase, ~] = fileparts(parts{idx+3});

% Assume the number is the first part of the file name, separated by an underscore.
fileParts = strsplit(fileBase, '_');
extractedNumber = fileParts{1};
end