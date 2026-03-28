%% Run one simulation: single candidate x single recording
% Returns the simulated c_cyto signal for one recording.
% This is the atomic unit of work for flattened parallelism.
function c_cyto = simulate_one_task(theta, other_settings, aqua_data_path, behavior, fps_org, fps_upsampled)
    run_seq      = behavior.run_seq;
    stim_onset   = behavior.stim_onset;
    reward_onset = behavior.reward_onset;
    TT           = behavior.TT;

    % precompute inputs (loads from cache if available)
    [pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun] = ...
        precompute_inputs(theta, other_settings, aqua_data_path, fps_org, fps_upsampled);

    % run simulation
    [chemical_table, ~] = single_round_simulation_one_recording(...
        theta, other_settings, run_seq, stim_onset, reward_onset, TT, ...
        pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun);
    c_cyto = chemical_table{:, "cytosolic_Ca"};
end
