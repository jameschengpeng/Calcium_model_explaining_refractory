function [X, names] = pop2mat(pop, selected_fields)
% POP2MAT   –  hyper‑parameter struct → numeric matrix
% Inputs:
%     pop   : 1xn   cell containers.Map
%     selected_fields: 1xdd cellstr. The fields that participate in
%     optimization
% Returns:
%     X     : n×d   numeric, n is the number of samples, d is
%     dimensionality of parameter space
%     names : 1×d   cellstr   (column order of X)

if nargin == 1
    selected_fields = pop{1}.keys;
end
n = numel(pop);
names = sort(selected_fields);     % fixed alphabetical order
d = numel(names); % dimensionality of selected-parameter space

X = zeros(n,d);
for i = 1:n
    for k = 1:d
        X(i,k) = pop{i}(names{k});
    end
end
end
