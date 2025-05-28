%  crossover: Combines two parent candidates to produce an offspring.
function offspring = crossover(parent1, parent2, important_fields)
offspring = copy_map(parent1);
% modify important_fields
for i = 1:length(important_fields)
    key = important_fields{i};
    % Uniform crossover: Randomly pick the parameter value from one of the parents.
    if rand < 0.5
        offspring(key) = parent1(key);
    else
        offspring(key) = parent2(key);
    end
end
end
