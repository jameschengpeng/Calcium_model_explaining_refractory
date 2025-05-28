%% precompute the input sequences and their functional form
function [precomputed_NE_seq, precomputed_NE_func, precomputed_glu_stim_seq, precomputed_glu_stim_func, precomputed_DA_seq, precomputed_DA_func] = ...
    precompute_inputs(theta, other_settings, aqua_data_path, fps_org, fps_upsampled)
% note that aqua_data_path is a mat file, remove the .mat will be a folder
recording_folder = aqua_data_path(1:end-4); % the folder to store the precomputed inputs
precomputed_input_file = fullfile(recording_folder, 'precomputed_inputs.mat');

if ~isfile(precomputed_input_file) % precomputed information is NOT available
    DA_spont = 0.01;
    running_threshold = 0.5;
    smoothing_factor = 5;
    NE_spont = 0.01;
    [~, run_seq, stim_onset, reward_onset, TT] = behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled);
    
    [precomputed_NE_seq, precomputed_NE_func] = run2NE_seq_func(run_seq, NE_spont, fps_org, fps_upsampled, TT, theta);

    [precomputed_glu_stim_seq, precomputed_glu_stim_func] = run2glu_stim_seq_func(run_seq, stim_onset, fps_org, fps_upsampled, theta('K_glu1'), theta('K_glu2'), other_settings, theta);
    
    [precomputed_DA_seq, precomputed_DA_func] = run2DA_seq_func(run_seq, DA_spont, reward_onset, fps_org, fps_upsampled, running_threshold, smoothing_factor, ...
        theta('K_DA_rew1'), theta('K_DA_rew2'), theta('b_DA_run'), theta('b_run'), theta('b_rew'));

    % save them for future use
    save(precomputed_input_file, "precomputed_NE_seq", "precomputed_NE_func", "precomputed_glu_stim_seq", "precomputed_glu_stim_func", "precomputed_DA_seq", "precomputed_DA_func");
else
    precomputed_inputs = load(precomputed_input_file);
    
    precomputed_NE_seq = precomputed_inputs.precomputed_NE_seq; precomputed_NE_func = precomputed_inputs.precomputed_NE_func;
    
    precomputed_glu_stim_seq = precomputed_inputs.precomputed_glu_stim_seq; precomputed_glu_stim_func = precomputed_inputs.precomputed_glu_stim_func;

    precomputed_DA_seq = precomputed_inputs.precomputed_DA_seq; precomputed_DA_func = precomputed_inputs.precomputed_DA_func;
end

end