function [rank, fronts] = nonDominatedSort(F, direction_indicator)
% nonDominatedSort performs fast non-dominated sorting on candidate solutions.
%
%   [rank, fronts] = nonDominatedSort(F)
%
% Inputs:
%   F - A matrix of objective values with dimensions m x n, where m is the
%       number of objectives and n is the number of candidates. Each column of F
%       is assumed to be a metric column vector for a candidate.
%   direction_indicator - A column vector of -1 and 1. Set 1 if this entry
%                         should be minimized, set -1 if this entry should
%                         be maximized
%
% Outputs:
%   rank   - A row vector of length n where rank(i) is the Pareto front number
%            (rank 1 is the best front) for the i-th candidate.
%   fronts - A cell array where fronts{j} contains the indices of candidates
%            belonging to the j-th Pareto front.
%
% Note: For this problem, for the first two metrics larger values are better,
%       and for the third metric smaller values are better.
%
F = round(F, 3); % precision up to 3 decimal places
nCandidates = size(F, 2);

% Initialize cell arrays and counters.
S = cell(nCandidates, 1);    % S{i} will contain indices of solutions dominated by candidate i.
n = zeros(1, nCandidates);   % n(i) counts the number of candidates that dominate candidate i.
rank = zeros(1, nCandidates);% The Pareto rank for each candidate.
fronts = {};                 % Cell array for grouping candidates by Pareto front.

% in case some old code just pass in F, the default metrics are corr,
% cosine and rmse. corr and cosine are to be maximized, while rmse is to be
% minimized. 
if nargin == 1
    if size(F, 1) == 3
        direction_indicator = [-1; -1; 1];
    elseif size(F, 1) == 4 % in this case, rmse for peaks is added, which should be minimized
        direction_indicator = [-1; -1; 1; 1];
    else
        error('Number of metrics not correct. It should be either 3 or 4');
    end
end


% Step 1: For each candidate, determine the set of dominated solutions and count dominations.
for p = 1:nCandidates
    S{p} = []; 
    n(p) = 0;
    for q = 1:nCandidates
        if p ~= q
            if dominates(F(:,p), F(:,q), direction_indicator)
                % Candidate p dominates candidate q.
                S{p} = [S{p}, q];
            elseif dominates(F(:,q), F(:,p), direction_indicator)
                % Candidate q dominates candidate p.
                n(p) = n(p) + 1;
            end
        end
    end
    % If candidate p is not dominated by any candidate, assign it to the first front.
    if n(p) == 0
        rank(p) = 1;
        if isempty(fronts) || length(fronts) < 1 || isempty(fronts{1})
            fronts{1} = p;
        else
            fronts{1} = [fronts{1}, p];
        end
    end
end

% Step 2: Build subsequent fronts.
i = 1;
while ~isempty(fronts{i})
    Q = []; % Temporary container for the next front.
    for p = fronts{i}
        for q = S{p}
            n(q) = n(q) - 1;
            if n(q) == 0
                rank(q) = i + 1;
                Q = [Q, q];
            end
        end
    end
    i = i + 1;
    fronts{i} = Q;
end

% Remove the last front if it is empty.
if isempty(fronts{end})
    fronts(end) = [];
end

end

function flag = dominates(a, b, direction_indicator)
% dominates checks if candidate a dominates candidate b.
%
%   flag = dominates(a, b)
%
% Both a and b are column vectors representing objective values.
% If an entry in direction_indicator is 1, then it should be minimized
% If an entry in direction_indicator is -1, then it should be maximized
%

flag = all(direction_indicator .* a <= direction_indicator .* b) && any(direction_indicator .* a < direction_indicator .* b);
end
