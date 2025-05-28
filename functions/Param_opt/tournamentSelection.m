%  Tournament Selection: Selects candidates from the population using tournament
%  selection based on Pareto rank.
function selected = tournamentSelection(population, ranks, tournament_size, num_selected)
pop_size = length(population);
selected = cell(num_selected, 1);
for i = 1:num_selected
    % Randomly choose tournament_size candidates.
    indices = randi(pop_size, tournament_size, 1);
    candidate_ranks = ranks(indices);
    [~, best_idx] = min(candidate_ranks);  % Lower rank is better.
    selected{i} = population{indices(best_idx)};
end
end