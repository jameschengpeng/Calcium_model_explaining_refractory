%% helper function
function [extractedName, extractedLocation, extractedNumber] = extractInfoFromPath(path)
    % Split path into parts using file separator
    parts = strsplit(path, filesep);

    % Extract required parts from the path
    extractedName = parts{end-3};      % 'Joey'
    extractedLocation = parts{end-2};  % 'RecLoc2'

    % Extract file name without extension
    [~, fileName, ~] = fileparts(path);

    % Extract number part (before underscore)
    underscoreIdx = strfind(fileName, '_');
    if ~isempty(underscoreIdx)
        extractedNumber = fileName(1:underscoreIdx(1)-1);
    else
        extractedNumber = fileName;
    end
end
