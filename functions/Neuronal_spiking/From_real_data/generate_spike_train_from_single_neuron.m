%% generate a spike train from spiking probability from one neuron
function spike_times = generate_spike_train_from_single_neuron(spike_prob_single, fps, prob_threshold)
%GENERATE_SPIKE_TRAIN_FROM_SINGLE_NEURON  Simulate spike times from per-bin probabilities
%   spike_prob_single : 1×N vector of P(spike in bin i), values in [0,1] or NaN
%   fps               : sampling rate in Hz (bins per second)
%   spike_times       : vector of spike times (in seconds)

% Impute missing probabilities (NaN) as 0
spike_prob_single(isnan(spike_prob_single)) = 0;

% keep only those above threshold
spike_prob_single(spike_prob_single < prob_threshold) = 0;

% Validate inputs
if any(spike_prob_single < 0) || any(spike_prob_single > 1)
    error('All entries of spike_prob_single must be between 0 and 1 (NaNs allowed).');
end
if ~isscalar(fps) || fps <= 0
    error('fps must be a positive scalar.');
end

% Bin width in seconds
dt = 1 / fps;

% Generate spikes by Bernoulli trials
isSpike = rand(size(spike_prob_single)) < spike_prob_single;

% Find the bin indices where spikes occurred
spike_bins = find(isSpike);

% Convert bin indices to times (seconds)
%   - If the first bin corresponds to time 0, spike at bin k occurs at (k-1)*dt
spike_times = (spike_bins - 1) * dt;
end
