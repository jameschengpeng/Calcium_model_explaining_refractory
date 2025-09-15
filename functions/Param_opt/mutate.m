%  mutate: Applies a random perturbation to a candidate's parameters.
function mutated_candidate = mutate(candidate, param_ranges, important_fields, mutation_rate)
% mutate with log-space awareness for positive parameters
% candidate: struct or containers.Map
% param_ranges.(key) = [lb, ub]

    mutated_candidate = candidate;

    % ----- choose one -----
    mode = 'log_uniform';  % 'log_uniform' or 'log_normal'
    log_step_decades = 0.2; % only used in 'log_normal' (~×10^0.2 ≈ ×1.58 per 1σ)

    for i = 1:length(important_fields)
        key = important_fields{i};
        if rand < mutation_rate
            range = param_ranges.(key);
            lb = range(1); ub = range(2);

            % Accessor depending on candidate type (struct vs. containers.Map)
            try
                x = candidate.(key);
                is_struct = true;
            catch
                x = candidate(key);
                is_struct = false;
            end

            if lb > 0 && ub > 0
                % --- LOG-SPACE MUTATION ---
                log_lb = log10(lb); log_ub = log10(ub);

                switch mode
                    case 'log_uniform'
                        % Global: resample uniformly in log10 space within [lb, ub]
                        log_x_new = log_lb + (log_ub - log_lb) * rand;
                    case 'log_normal'
                        % Local: multiplicative Gaussian step in log10 space
                        log_x_cur = log10(max(min(x, ub), lb));
                        log_x_new = log_x_cur + log_step_decades * randn;
                        % keep inside bounds
                        log_x_new = min(max(log_x_new, log_lb), log_ub);
                    otherwise
                        error('Unknown mode');
                end

                mutated_value = 10.^log_x_new;

            else
                % --- FALLBACK: LINEAR ADDITIVE MUTATION (original style) ---
                perturbation = (ub - lb) * 0.1 * (rand - 0.5);
                mutated_value = x + perturbation;
            end

            % Clamp to bounds
            mutated_value = max(lb, min(ub, mutated_value));

            % Write back
            if is_struct
                mutated_candidate.(key) = mutated_value;
            else
                mutated_candidate(key) = mutated_value;
            end
        end
    end
end


% function mutated_candidate = mutate(candidate, param_ranges, important_fields, mutation_rate)
% mutated_candidate = candidate;
% for i = 1:length(important_fields)
%     key = important_fields{i};
%     if rand < mutation_rate
%         range = param_ranges.(key);
%         % Perturbation scaled as 10% of the parameter's range.
%         perturbation = (range(2)-range(1)) * 0.1 * (rand - 0.5);
%         mutated_value = candidate(key) + perturbation;
%         % Ensure the mutated value is within bounds.
%         mutated_candidate(key) = max(range(1), min(range(2), mutated_value));
%     end
% end
% end