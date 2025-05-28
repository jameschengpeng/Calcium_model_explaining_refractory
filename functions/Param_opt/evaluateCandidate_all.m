%% to evaluate a candidate parameter configuration on all recordings
function metrics = evaluateCandidate_all(candidate, other_settings, data_path, fps_org, fps_upsampled)
    % Run simulation; inner parfor is active within the simulation function.
    save_tables = false;
    important_fields = false;
    [sim_result, real_result] = single_round_simulation_all_recordings(candidate, other_settings, data_path, fps_org, fps_upsampled, save_tables, important_fields);
    % Compute multi-objective metrics.
    metrics = multi_metric_eval(sim_result, real_result);
end