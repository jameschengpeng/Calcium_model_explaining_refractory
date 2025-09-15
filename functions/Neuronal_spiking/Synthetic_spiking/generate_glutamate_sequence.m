function [glu_trace, u_trace, x_trace, spike_times] = generate_glutamate_sequence(rate_hz, n_steps, fps, U, tau_F, tau_D, tau_glu)
duration_sec = 1/fps * (n_steps-1);
spike_times = generate_spike_train(rate_hz, duration_sec);
[glu_trace, u_trace, x_trace] = spike_train_to_glu(spike_times, n_steps, fps, U, tau_F, tau_D, tau_glu);
end
