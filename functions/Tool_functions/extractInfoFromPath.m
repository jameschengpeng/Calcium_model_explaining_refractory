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