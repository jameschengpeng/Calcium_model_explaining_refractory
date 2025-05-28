%% preprocessing of the data
function [resp_seq, run_seq, stim_onset, reward_onset, TT, regional_response_temp_down] = behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled)
aqua_data = load(aqua_data_path);
TT = aqua_data.TT;
full_vel = aqua_data.velocity; % full length velocity
regional_response_temp_down = aqua_data.avg_roi_temp_down; % temporally downsampled

[regional_response_same_baseline, ~] = response_processing(regional_response_temp_down, 3, 100);
regional_response_same_baseline = series_augmentation(regional_response_same_baseline, full_vel);

run_seq = abs(full_vel);
resp_seq = scale_resp(regional_response_same_baseline); % rescale the response sequence
resp_func = get_resp_func(resp_seq, fps_org); % resp_func = @(t)
termination_time = (length(resp_seq)-1)/fps_org;
n_steps = ceil(fps_upsampled * termination_time);
discrete_time = linspace(0, termination_time, n_steps);
resp_seq = arrayfun(resp_func, discrete_time);

visual_run_stim_info = TT{:, "visual_run_stim"};
visual_sit_stim_info = TT{:, "visual_sit_stim"};
trial_type = TT{:, "trial_type"};
visual_stim_info = visual_run_stim_info; % only the stim for running counts for visual evoked glutamate
for ii = 1:length(visual_run_stim_info)
    if trial_type(ii) == 3 % correct rejection
        visual_stim_info(ii) = visual_sit_stim_info(ii);
    end
end

reward_info = TT{:, "reward"};
stim_onset = [];
for ii = 1:length(visual_stim_info)
    if visual_stim_info(ii) ~= -1 && TT{ii, "salient_stim"} ~= -1
        stim_onset = [stim_onset; visual_stim_info(ii)];
    elseif visual_stim_info(ii) ~= -1 && TT{ii, "threshold_stim"} ~= -1
        stim_onset = [stim_onset; -visual_stim_info(ii)];
    elseif visual_stim_info(ii) ~= -1 && TT{ii, "trial_type"} == 3 % for CR specially
        stim_onset = [stim_onset; visual_stim_info(ii)];
    end
end

reward_onset = reward_info(reward_info ~= -1);



end