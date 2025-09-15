% transforms the spike train to glutamate time series
function [glu_trace, u_trace, x_trace] = spike_train_to_glu(spike_times, n_steps, fps, U, tau_F, tau_D, tau_glu)
dt = 1/fps;
glu_trace = zeros(1, n_steps);
u_trace = zeros(1, n_steps); u_trace(1) = U; % utilization of synapses
x_trace = zeros(1, n_steps); x_trace(1) = 1; % percentage of available resources

spike_idx = 1;
n_spikes = length(spike_times);
t_spike = spike_times(spike_idx);

u_curr = u_trace(1);
x_curr = x_trace(1);
glu_curr = glu_trace(1);

% Define the derivative functions for the continuous part (except dirac
% delta)
f_u = @(u) -(u - U)/tau_F;
f_x = @(x) (1 - x)/tau_D;
f_glu = @(glu) -glu/tau_glu;

for ii = 2:n_steps
    t_curr = (ii-1) * dt;
    t_prev = (ii-2) * dt;
    if (spike_idx <= n_spikes) && (t_spike > t_prev) && (t_spike <= t_curr) % there is AT LEAST ONE spike in this interval
        while spike_idx <= n_spikes && t_spike <= t_curr % iterate all spikes within this interval
            r = u_curr * x_curr;
            % before the RK4 step, apply the spike jump at time t_spike
            u_curr = u_curr + U * (1-u_curr);
            x_curr = x_curr - u_curr * x_curr;
            glu_curr = glu_curr + r;
    
            spike_idx = spike_idx + 1;
            if spike_idx <= n_spikes
                t_spike = spike_times(spike_idx);
            end
        end
    end

    % RK4 updates
    % For u:
    k1_u = dt * f_u(u_curr);
    k2_u = dt * f_u(u_curr + 0.5*k1_u);
    k3_u = dt * f_u(u_curr + 0.5*k2_u);
    k4_u = dt * f_u(u_curr + k3_u);
    u_next = u_curr + (k1_u + 2*k2_u + 2*k3_u + k4_u) / 6;
    
    % For x:
    k1_x = dt * f_x(x_curr);
    k2_x = dt * f_x(x_curr + 0.5*k1_x);
    k3_x = dt * f_x(x_curr + 0.5*k2_x);
    k4_x = dt * f_x(x_curr + k3_x);
    x_next = x_curr + (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6;
    
    % For glu:
    k1_glu = dt * f_glu(glu_curr);
    k2_glu = dt * f_glu(glu_curr + 0.5*k1_glu);
    k3_glu = dt * f_glu(glu_curr + 0.5*k2_glu);
    k4_glu = dt * f_glu(glu_curr + k3_glu);
    glu_next = glu_curr + (k1_glu + 2*k2_glu + 2*k3_glu + k4_glu) / 6;

    % Save new values into arrays:
    u_trace(ii) = u_next;
    x_trace(ii) = x_next;
    glu_trace(ii) = glu_next;
    
    % Update current state for the next step:
    u_curr = u_next;
    x_curr = x_next;
    glu_curr = glu_next;
end

end


