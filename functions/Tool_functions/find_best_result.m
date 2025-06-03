% for manual tuned results, retrieve the one with largest correlation
% coefficient
function [best_result_filepath, best_result_corr] = find_best_result(path, prefix)
    matFiles = dir(fullfile(path, '*.mat'));
    fileNames = {matFiles.name};
    correlations = zeros(1, length(fileNames));

    for ii = 1:length(fileNames)
        file = fileNames{ii};

        % Check if the file name starts with "Results"
        if startsWith(file, prefix)
            decimalNumber = regexp(file, strcat(prefix, '_([\d\.]+)\.mat'), 'tokens');
            if ~isempty(decimalNumber)
                corr = str2double(decimalNumber{1}{1});
                correlations(ii) = corr;
            else
                correlations(ii) = -Inf; % If the file doesn't match the pattern, set correlation to -Inf
            end
        else
            correlations(ii) = -Inf; % Skip files that don't start with "Results"
        end
    end

    [best_result_corr, max_corr_idx] = max(correlations);
    best_result_filepath = fullfile(path, fileNames{max_corr_idx});
end