%% retrieve the search result of the genetic algorithm
% find_last can either be boolean, or a specific number of generation
function [result, idx] = find_GA_search_result(save_path, find_last)
% List matching files
pattern = fullfile(save_path, 'Generation_*.mat');
d = dir(pattern);
if isempty(d)
    error('No files matching ''Generation_*.mat'' found in %s.', save_path);
end

names = {d.name};
nums  = zeros(size(names));

% Extract the number after the underscore, before ".mat"
for k = 1:numel(names)
    tok = regexp(names{k}, '^Generation_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tok)
        error('Filename "%s" does not match ''Generation_xx.mat''.', names{k});
    end
    nums(k) = str2double(tok{1});
end

% Pick the first or last generation
if islogical(find_last) % extract either gen 0 or gen max
    if find_last
        idx = max(nums);
    else
        idx = min(nums);
    end
else % extract the specific gen
    idx = find_last;
end
filename = strcat('Generation_', num2str(idx), '.mat');

% Return full path
result = fullfile(save_path, filename);
end
