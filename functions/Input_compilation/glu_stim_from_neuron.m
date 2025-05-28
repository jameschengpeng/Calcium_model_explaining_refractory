%% generate the glu_stim sequence and function from neuronal spiking probability
function [glu_stim_seq, glu_stim_func] = glu_stim_from_neuron(spike_prob_visual, spike_prob_run, fps_org, fps_upsampled)
n_neuron = size(spike_prob_visual, 1) + size(spike_prob_run, 1);
n_visual = size(spike_prob_visual, 1);
n_steps = size(spike_prob_visual, 2);
glu_stim_seq_multi = zeros(n_neuron, n_steps);

U = 0.2;
tau_F = 1.5;
tau_D = 0.2;
tau_glu = 0.5;

for ii = 1:n_neuron
    if ii <= n_visual
        spike_prob_single = spike_prob_visual(ii, :);
    else
        idx = ii-n_visual;
        spike_prob_single = spike_prob_run(idx, :);
    end
    prob_threshold = 0.1;
    spike_times = generate_spike_train_from_single_neuron(spike_prob_single, fps_org, prob_threshold);
    [glu_trace, ~, ~] = spike_train_to_glu(spike_times, n_steps, fps_org, U, tau_F, tau_D, tau_glu);
    glu_stim_seq_multi(ii, :) = glu_trace;
end
% now, still in fps_org
% glu_stim_seq = mean(glu_stim_seq_multi, 1);
w_visual = 0.9;
avg_visual = mean(glu_stim_seq_multi(1:n_visual, :), 1);   % 1×N
avg_run = mean(glu_stim_seq_multi(n_visual+1:end, :), 1);   % 1×N
% weighted combination (still a 1×N row)
glu_stim_seq = w_visual * avg_visual + (1 - w_visual) * avg_run;
glu_stim_seq = scale_resp(glu_stim_seq, 0, 1);

termination_time = (n_steps-1)/fps_org;
discrete_time = linspace(0, termination_time, n_steps);
glu_stim_func = @(t) interp1(discrete_time, glu_stim_seq, t, 'linear', 'extrap');

% finally, upsample glu_stim_seq
n_steps_upsampled = ceil(fps_upsampled * termination_time);
discrete_time_upsampled = linspace(0, termination_time, n_steps_upsampled);
glu_stim_seq = arrayfun(glu_stim_func, discrete_time_upsampled);

end