%% given the run_seq & visual stimuli, compile a sequence & function of stimulation which leads to neuronal glutamate synthesis
% at time (s) stim_onset, visual stimulation is presented
% similar to NE, make the glutamate level in range [0,1]
function [glu_stim_seq, glu_stim_func] = run2glu_stim_seq_func(run_seq, stim_onset, fps_org, fps_upsampled, K_glu1, K_glu2, other_settings, theta)
NT_delay = theta('NT_delay');
% preprocess the run_seq in case fps_upsampled > 30glu_seq
run_func = seq2run_func(run_seq, fps_org); % run_func = @(t)
termination_time = (length(run_seq)-1)/fps_org;
n_steps = ceil(fps_upsampled * termination_time);
discrete_time = linspace(0, termination_time, n_steps);
run_seq = arrayfun(run_func, discrete_time); % now the run_seq is scaled

glu_max1 = 1.2; % for salient stim
glu_max2 = 1; % for threshold stim
%% visual stimuli evoked
visual_glu_stim = zeros(size(run_seq));
if ~isempty(stim_onset)
    % model the glutamate level after visual stimulation by product of two
    % Hill equations
    for ii = 1:length(stim_onset)
        if stim_onset(ii) > 0 % salient stim
            onset_idx = stim_onset(ii);
            glu_max = glu_max1;
        else % threshold stim
            onset_idx = -stim_onset(ii);
            glu_max = glu_max2;
        end
        stim_idx = round(onset_idx * fps_upsampled/fps_org); % rescale the stim_idx to the augmented domain
        local_length = min(fps_upsampled*20, length(visual_glu_stim)-stim_idx+1);
        local_glu_stim_seq = zeros(1, local_length); % visual stim induced glutamate lapse for some time
        glu_amp = glu_max;
        func = @(t) Hill_func(max(0, t-NT_delay), 1, K_glu1) * reverse_Hill_func(max(0, t-NT_delay), 1, K_glu2);
        t_seq = (0:(length(local_glu_stim_seq)-1))/fps_upsampled;
        excess_glu = arrayfun(func, t_seq);
        excess_glu = (excess_glu - min(excess_glu))/(max(excess_glu) - min(excess_glu)) * glu_amp;
        visual_glu_stim(stim_idx:(stim_idx+length(local_glu_stim_seq)-1)) = visual_glu_stim(stim_idx:(stim_idx+length(local_glu_stim_seq)-1)) + excess_glu;
    end
    visual_glu_stim = (visual_glu_stim - min(visual_glu_stim))./(max(visual_glu_stim) - min(visual_glu_stim)); % normalize to range 0~1
end

%% running evoked
run_seq_filt = imgaussfilt(run_seq, fps_upsampled); % apply gaussian filtering to obtain lower freq components
max_velocity = max(run_seq_filt);
running_periods = run_seq_filt >= 0.5;
start_indices = find(diff([0, running_periods]) == 1); % start points of running period
end_indices = find(diff([running_periods, 0]) == -1); % end points of running period
% only keep the periods >= fps_upsampled steps
revised_start = []; revised_end = [];
for ii = 1:length(start_indices)
    if end_indices(ii) - start_indices(ii) > fps_upsampled
        revised_start = [revised_start start_indices(ii)];
        revised_end = [revised_end end_indices(ii)];
    end
end
start_indices = revised_start; end_indices = revised_end;
run_glu_stim = zeros(size(run_seq));

for ii = 1:length(start_indices)
    s = start_indices(ii);
    e = end_indices(ii);
    time_length = (e-s)/fps_upsampled;
    run_ampl = max(run_seq(s:e));
    glu_stim_Hill = @(x) Hill_func(x, 2, 0.1*time_length) * reverse_Hill_func(x, 2, 0.3*time_length);
    part1 = zeros(1,s); part2 = (1:(length(run_seq)-s))./fps_upsampled;
    glu_stim_subseq = [part1 part2];
    glu_stim_subseq = arrayfun(glu_stim_Hill, glu_stim_subseq);
    max_val = max(glu_stim_subseq);
    glu_stim_subseq = glu_stim_subseq./max_val;
    glu_stim_subseq = glu_stim_subseq.*(run_ampl/max_velocity);
    run_glu_stim = run_glu_stim + glu_stim_subseq;
end

run_glu_stim = (run_glu_stim - min(run_glu_stim))./(max(run_glu_stim) - min(run_glu_stim));

%% combining visual evoked and run evoked glu_stim
glu_stim_seq = 0.8 * run_glu_stim + 0.2 * visual_glu_stim;
% glu_stim_seq = 0.3 * run_glu_stim + 0.7 * visual_glu_stim;
% for glu_stim_func, which should take continuous input t, do linear
% interpolation to the glu_stim_seq
time_course = (0:(length(glu_stim_seq)-1))/fps_upsampled;

% glu_stim_func = @(t) interp1(time_course, glu_stim_seq, t, 'linear', 'extrap');
F = griddedInterpolant(time_course, glu_stim_seq, 'linear', 'linear');
glu_stim_func = @(t) F(t);

% when simulating the case without the participation of glutamate
if ~other_settings('glutamate')
    glu_stim_seq = zeros(size(glu_stim_seq));
    glu_stim_func = @(t) interp1(time_course, glu_stim_seq, t, 'linear', 'extrap');
end
end