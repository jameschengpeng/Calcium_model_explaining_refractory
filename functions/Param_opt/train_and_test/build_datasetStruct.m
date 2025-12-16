%% to create dataset struct, key is mouse's name, value is a cell array containing the path to the recording 
% Input: data_path   the path to the Data folder
% Output: datasetStruct   to be fed to the function
% createStratifiedTrainTestSplits()
function datasetStruct = build_datasetStruct(data_path)
mice_names = {'Joey', 'Ross', 'Pecan', 'Peanut'};
% build a list of dataset struct
datasetStruct = struct();

for ii = 1:length(mice_names)
    mouse = mice_names{ii};
    datasetStruct.(mouse) = {}; % use a cell array to store the filepath for recordings to each mouse
    RecLocs = find_subDirs(fullfile(data_path, mouse));
    for jj = 1:length(RecLocs)
        location = RecLocs{jj};
        location_path = fullfile(data_path, mouse, location);
        recordings = find_subDirs(location_path);
        unique_prefix = find_unique_recording_prefix(recordings);
        for kk = 1:length(unique_prefix)
            prefix = unique_prefix{kk};
            
            dirPattern = fullfile(location_path, [prefix, '*']);

            % List the matching directory
            matched_dir = dir(dirPattern);
            
            % the matFile has the same name as the matched_dir
            matFile = strcat(matched_dir.name, '.mat');
            matFilePath = fullfile(matched_dir.folder, matched_dir.name, matFile);

            % Construct the full path to the selected file.
            datasetStruct.(mouse) = [datasetStruct.(mouse); matFilePath];
        end
    end
end

end

function unique_recording_prefix = find_unique_recording_prefix(recordings)
prefixes = cellfun(@(s) s(1:strfind(s, '_')-1), recordings, 'UniformOutput', false);
unique_recording_prefix = unique(prefixes);
end