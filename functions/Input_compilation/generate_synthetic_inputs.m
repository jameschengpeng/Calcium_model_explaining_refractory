%% generate synthetic run_seq, stim_onset, and reward_onset. sim for simulated
% NE_time, glu_time, DA_time should all be in seconds
% all of them are 2*n, 1st row for starting time, 2nd row for ending time
function [sim_NE_seq, sim_NE_func, sim_glu_stim_seq, sim_glu_stim_func, sim_DA_seq, sim_DA_func] = ...
    generate_synthetic_inputs(fps_upsampled, time_length_seconds, NE_time, glu_time, DA_time, NE_spont, DA_spont)
n_steps = ceil(time_length_seconds * fps_upsampled);
sim_NE_seq = NE_spont .* ones(1, n_steps);
sim_glu_stim_seq = 0.01 * ones(1, n_steps);
sim_DA_seq = DA_spont .* ones(1, n_steps);

sim_NE_seq = simulate_signal_series(NE_time, fps_upsampled, sim_NE_seq, false);
% sim_glu_stim_seq = simulate_signal_series_from_spikes(glu_time, fps_upsampled, sim_glu_stim_seq);
sim_glu_stim_seq = simulate_signal_series(glu_time, fps_upsampled, sim_glu_stim_seq, true);
sim_DA_seq = simulate_signal_series(DA_time, fps_upsampled, sim_DA_seq, true);

sim_NE_func = seq2func(sim_NE_seq, fps_upsampled);
sim_glu_stim_func = seq2func(sim_glu_stim_seq, fps_upsampled);
sim_DA_func = seq2func(sim_DA_seq, fps_upsampled);
end

%% tool functions
% generate simulated glutamate trace from neuronal spiking
function sim_seq = simulate_signal_series_from_spikes(NT_time, fps, sim_seq)
rate_hz = 5;
U = 0.2;
tau_F = 1.5;
tau_D = 0.2;
tau_glu = 0.5;
if ~ isempty(NT_time)
    for ii = 1:size(NT_time, 2)
        start_time = NT_time(1, ii);
        end_time = NT_time(2, ii);
        n_steps = (end_time-start_time) * fps;
        [signal_part,~,~,~] = generate_glutamate_sequence(rate_hz, n_steps, fps, U, tau_F, tau_D, tau_glu);
        start_idx = ceil(start_time * fps);
        sim_seq(start_idx:(start_idx+length(signal_part)-1)) = sim_seq(start_idx:(start_idx+length(signal_part)-1)) + signal_part;
    end
end
sim_seq = scale_resp(sim_seq, 0.01, 1.2);
end

function sim_seq = simulate_signal_series(NT_time, fps, sim_seq, skew)
if ~ isempty(NT_time)
    for ii = 1:size(NT_time, 2)
        start_time = NT_time(1, ii);
        end_time = NT_time(2, ii);
        if ~ skew
            signal_peak = simulate_signal_series_one_peak_no_skew(start_time, end_time, fps);
        else
            signal_peak = simulate_signal_series_one_peak_skew(start_time, end_time, fps);
        end
        start_idx = ceil(start_time * fps);
        sim_seq(start_idx:(start_idx+length(signal_peak)-1)) = sim_seq(start_idx:(start_idx+length(signal_peak)-1)) + signal_peak;
    end
end
end

function signal = simulate_signal_series_one_peak_no_skew(start_time, end_time, fps)
t = start_time:1/fps:end_time;
mu = (start_time + end_time) / 2;
sigma = (end_time - start_time) / 5;
signal = exp(-(t-mu).^2./(2*sigma.^2));
end

function signal = simulate_signal_series_one_peak_skew(start_time, end_time, fps)
t = start_time:1/fps:end_time;
mu = (start_time + end_time) / 2;
sigma = (end_time - start_time) / 5;
signal = exp(-(t-mu).^2./(2*sigma.^2)) .* (1 + erf(10*(t-mu)./(sqrt(2).*sigma))) ./ 2;
end

function seq_func = seq2func(seq, fps)
time_course = (0:(length(seq)-1))/fps;
seq_func = @(t) interp1(time_course, seq, t, 'linear', 'extrap');
end