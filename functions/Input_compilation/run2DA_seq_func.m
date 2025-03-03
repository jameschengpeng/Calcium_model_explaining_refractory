%% given the run_seq, reward, recover the DA_seq and DA_func
% similar to NE, make DA at level [0,1]
function [DA_seq, DA_func] = run2DA_seq_func(run_seq, DA_spont, reward_onset, fps_org, fps_upsampled, running_threshold, smoothing_factor, ...
    K_DA_rew1, K_DA_rew2, b_DA_run, b_run, b_rew)
% preprocess run_seq in case fps_upsampled > fps_org
run_func = seq2run_func(run_seq, fps_org); % run_func = @(t)
termination_time = (length(run_seq)-1)/fps_org;
n_steps = ceil(fps_upsampled * termination_time);
discrete_time = linspace(0, termination_time, n_steps);
run_seq = arrayfun(run_func, discrete_time);

%% produce the running DA
run_seq_filt = imgaussfilt(run_seq, smoothing_factor); % apply gaussian filtering to obtain lower freq components
max_velocity = max(run_seq_filt);
running_periods = run_seq_filt >= running_threshold;
start_indices = find(diff([0, running_periods]) == 1); % start points of running period
end_indices = find(diff([running_periods, 0]) == -1); % end points of running period
% only keep the periods >= 3 * fps_upsampled steps, i.e., 1s
revised_start = []; revised_end = [];
for ii = 1:length(start_indices)
    if end_indices(ii) - start_indices(ii) > 3 * fps_upsampled
        revised_start = [revised_start start_indices(ii)];
        revised_end = [revised_end end_indices(ii)];
    end
end
start_indices = revised_start; end_indices = revised_end;
run_DA = zeros(size(run_seq));

for ii = 1:length(start_indices)
    s = start_indices(ii);
    e = end_indices(ii);
    time_length = (e-s)/fps_upsampled;
    run_ampl = max(run_seq(s:e));
    DA_Hill = @(x) Hill_func(x, 2, 0.1*time_length) * reverse_Hill_func(x, 2, 0.3*time_length);
    part1 = zeros(1,s); part2 = (1:(length(run_seq)-s))./fps_upsampled;
    DA_seq = [part1 part2];
    DA_seq = arrayfun(DA_Hill, DA_seq);
    max_val = max(DA_seq);
    DA_seq = (DA_seq./max_val).*(1-DA_spont);
    DA_seq = DA_seq.*(run_ampl/max_velocity);
    run_DA = run_DA + DA_seq;
end

run_DA = (run_DA - min(run_DA))/(max(run_DA) - min(run_DA));

%% produce the reward DA
reward_DA = DA_spont * ones(size(run_seq));
if isempty(reward_onset)
    DA_seq = b_rew * reward_DA + b_run * run_DA;
else
    % model the dopamine level after given the reward by product of two
    % Hill equations
    for ii = 1:length(reward_onset)
        reward_idx = round(reward_onset(ii)/fps_org * fps_upsampled); % rescale the reward_idx to the augmented domain
        DA_subseq = DA_spont * ones(1, (fps_upsampled*20));
        c1 = K_DA_rew1; c2 = K_DA_rew2;
        func = @(t) t^2/(t^2+c1^2) * c2^2/(c2^2 + t^2);
        t_seq = 1:(length(DA_subseq));
        t_seq = t_seq./fps_upsampled;
        end_idx = min(length(run_seq), reward_idx+length(DA_subseq)-1);
        len_reward_DA = end_idx - reward_idx + 1;
        reward_DA(reward_idx:end_idx) = arrayfun(func, t_seq(1:len_reward_DA));
    end
    reward_DA = (reward_DA - min(reward_DA))/(max(reward_DA) - min(reward_DA));
    DA_seq = b_rew * reward_DA + b_run * run_DA;
end

if isempty(reward_onset)
    DA_seq = b_run * (DA_seq - min(DA_seq))/(max(DA_seq) - min(DA_seq));
else
    DA_seq = (DA_seq - min(DA_seq))/(max(DA_seq) - min(DA_seq));
end

% for DA_func, which should take continuous input t, do linear
% interpolation to the DA_seq
time_course = (1:length(DA_seq))/fps_upsampled;
DA_func = @(t) interp1(time_course, DA_seq, t, 'linear', 'extrap');
end