% given the spiking rate (in Hz) and the total duration in seconds, return
% an array of time points for neuronal spiking
function spike_times = generate_spike_train(rate_hz, duration_sec)
% rate_hz: firing rate (spikes/sec)
% duration_sec: total duration of spike train (sec)

spike_times = [];
t = 0;
while t < duration_sec
    isi = exprnd(1/rate_hz);  % exponentially-distributed inter-spike intervals
    t = t + isi;
    if t <= duration_sec
        spike_times = [spike_times; t];
    end
end
end
