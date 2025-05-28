%% to evaluate a candidate parameter configuration on subset of recordings
function metrics = evaluateCandidate_part(candidate, other_settings, train_or_test_set, fps_org, fps_upsampled, save_path)
    % Run simulation; inner parfor is active within the simulation function.
    save_tables = false;
    important_fields = false;
    [sim_result, real_result] = single_round_simulation_part_recordings(...
        candidate, other_settings, train_or_test_set, fps_org, fps_upsampled, save_tables, important_fields);
    % Compute multi-objective metrics.
    metrics = multi_metric_eval(sim_result, real_result, save_path);
end