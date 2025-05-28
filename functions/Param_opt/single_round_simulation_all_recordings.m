%% for a given theta, simulate on all recordings
function [simulation_result, real_result] = single_round_simulation_all_recordings(...
    theta, other_settings, data_path, fps_org, fps_upsampled, save_tables, important_fields)

% — get all the .mat paths —
aqua_data_path_ALL = get_aqua_data_path_ALL(data_path);
numFiles = numel(aqua_data_path_ALL);

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
    aqua_data_path = aqua_data_path_ALL{ii};
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
            fname = sprintf('GA_Results_important_fields_%.3f.mat', corr_val);
        else
            fname = sprintf('GA_Results_all_fields_%.3f.mat', corr_val);
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
    path = aqua_data_path_ALL{ii};
    [nm, loc, num] = extractInfoFromPath(path);
    field = sprintf('%s_%s_%s', nm, loc, num);
    simulation_result.(field) = simulation_result_tmp{ii};
    real_result.(field)       = real_result_tmp{ii};
end

end


%% helper functions
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



function aqua_data_path_ALL = get_aqua_data_path_ALL(data_path)
mice_names = {'Joey', 'Ross', 'Pecan', 'Peanut'};
% build a list of aqua_data_path for the recordings
aqua_data_path_ALL = {};
for ii = 1:length(mice_names)
    mouse = mice_names{ii};
    RecLocs = find_subDirs(fullfile(data_path, mouse));
    for jj = 1:length(RecLocs)
        location = RecLocs{jj};
        location_path = fullfile(data_path, mouse, location);
        recordings = find_subDirs(location_path);
        unique_prefix = find_unique_recording_prefix(recordings);
        for kk = 1:length(unique_prefix)
            prefix = unique_prefix{kk};
            filePattern = fullfile(location_path, [prefix, '*.mat']);

            % List all the matching files
            allFiles = dir(filePattern);
            
            % Exclude files that have "processed" in their filename (case insensitive)
            mat_files = allFiles(~contains({allFiles.name}, 'processed', 'IgnoreCase', true));
            
            % Check if any file was found
            if isempty(mat_files)
                error('No MAT file found with the specified prefix.');
            end
            
            % Extract the names from the file structure into a cell array.
            fileNames = {mat_files.name};
            
            % Sort the file names lexicographically (i.e., based on ASCII code order).
            [~, idx] = sort(fileNames);
            
            % The first element in the sorted list is the file with the smallest ASCII code
            selectedFile = mat_files(idx(1));
            
            % Construct the full path to the selected file.
            selectedFilePath = fullfile(selectedFile.folder, selectedFile.name);
            aqua_data_path_ALL{end+1} = selectedFilePath;
        end
    end
end
end

function unique_recording_prefix = find_unique_recording_prefix(recordings)
prefixes = cellfun(@(s) s(1:strfind(s, '_')-1), recordings, 'UniformOutput', false);
unique_recording_prefix = unique(prefixes);
end

function subDirs = find_subDirs(path)
allFiles = dir(path);
isSubDir = [allFiles.isdir];
subDirs = {allFiles(isSubDir).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'}));
end

function parsaveStruct(filename, S)
    % S must now be a 1×1 struct
    assert( isstruct(S) && numel(S)==1, ...
        'parsaveStruct: need a scalar struct' );

    % '-struct' (or alias '-fromstruct') takes the fields of S
    % and writes them as variables into the .mat
    save(filename, '-struct', 'S', '-v7.3')
end
