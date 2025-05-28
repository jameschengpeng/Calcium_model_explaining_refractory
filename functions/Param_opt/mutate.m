%  mutate: Applies a random perturbation to a candidate's parameters.
function mutated_candidate = mutate(candidate, param_ranges, important_fields, mutation_rate)
mutated_candidate = candidate;
for i = 1:length(important_fields)
    key = important_fields{i};
    if rand < mutation_rate
        range = param_ranges.(key);
        % Perturbation scaled as 10% of the parameter's range.
        perturbation = (range(2)-range(1)) * 0.1 * (rand - 0.5);
        mutated_value = candidate(key) + perturbation;
        % Ensure the mutated value is within bounds.
        mutated_candidate(key) = max(range(1), min(range(2), mutated_value));
    end
end
end