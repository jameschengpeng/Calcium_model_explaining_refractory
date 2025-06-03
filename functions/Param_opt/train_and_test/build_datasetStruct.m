%% to create dataset struct, key is mouse's name, value is a cell array containing the path to the recording 
% Input: data_path   the path to the folder of AQUA_processed_data
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
            filePattern = fullfile(location_path, [prefix, '*'], [prefix, '*.mat']);

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
            datasetStruct.(mouse) = [datasetStruct.(mouse); selectedFilePath];
        end
    end
end

end

function unique_recording_prefix = find_unique_recording_prefix(recordings)
prefixes = cellfun(@(s) s(1:strfind(s, '_')-1), recordings, 'UniformOutput', false);
unique_recording_prefix = unique(prefixes);
end