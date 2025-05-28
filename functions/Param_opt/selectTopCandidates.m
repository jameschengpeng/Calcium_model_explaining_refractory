%  selectTopCandidates: Selects the best pop_size candidates from the merged
%  population based on Pareto ranks and also returns their fitness values.
function [selected, selected_fitness] = selectTopCandidates(merged_population, merged_ranks, F_merged, pop_size)
% Sort indices by ascending Pareto rank.
[~, sort_idx] = sort(merged_ranks);
selected = merged_population(sort_idx(1:pop_size));
selected_fitness = F_merged(:, sort_idx(1:pop_size));
end