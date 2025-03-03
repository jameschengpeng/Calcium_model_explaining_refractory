%% given the run_seq, recover the NE_seq and the NE_func
function [NE_seq, NE_func] = run2NE_seq_func(run_seq, NE_spont, fps_org, fps_upsampled, TT, theta)
NT_delay = theta('NT_delay');
tau = theta('leaky_factor');
boost_factor = 2;
run_func = seq2run_func(run_seq, fps_org); % run_func = @(t)
termination_time = (length(run_seq)-1)/fps_org;
n_steps = ceil(fps_upsampled * termination_time);
discrete_time = linspace(0, termination_time, n_steps);
run_seq = arrayfun(run_func, discrete_time);
NE_seq = NE_spont * ones(1, length(run_seq));
booster = CR_booster(boost_factor, TT, fps_org, fps_upsampled, NE_seq);
dt = 1/fps_upsampled;
NE = NE_spont;
for ii = (0*fps_upsampled):(length(run_seq)-1)
    t = max(0, ii/fps_upsampled - NT_delay);
    boost_coeff = booster(ii+1);
    dNE_dt_func = get_dNE_dt_func();
    k1 = dt * dNE_dt_func(NE, run_func(t), tau, boost_coeff);
    k2 = dt * dNE_dt_func(NE+1/2*k1, run_func(t+dt/2), tau, boost_coeff);
    k3 = dt * dNE_dt_func(NE+1/2*k2, run_func(t+dt/2), tau, boost_coeff);
    k4 = dt * dNE_dt_func(NE+k3, run_func(t+dt), tau, boost_coeff);
    NE = NE + k1/6 + k2/3 + k3/3 + k4/6;
    NE_seq(ii+1) = NE;
end
NE_seq = NE_seq./max(NE_seq);
NE_seq = NE_seq.^1;
time_course = (0:(length(NE_seq)-1))/fps_upsampled;

% NE_func = @(t) interp1(time_course, NE_seq, t, 'linear', 'extrap');

% Create the interpolant outside of the function handle
F = griddedInterpolant(time_course, NE_seq, 'linear', 'linear');  % 'linear' for extrapolation mode too
NE_func = @(t) F(t);
end

function dNE_dt_func = get_dNE_dt_func()
dNE_dt_func = @(NE, vel, tau, boost_coeff) 1/tau * (-NE + boost_coeff * vel);
end

function booster = CR_booster(boost_factor, TT, fps_org, fps_upsampled, NE_seq)
booster = ones(size(NE_seq));
visual_sit_stim = TT{:, "visual_sit_stim"};
trial_type = TT{:, "trial_type"};
CR_stim = visual_sit_stim(trial_type == 3);
CR_ind = ceil(CR_stim .* (fps_upsampled/fps_org));
for ii = 1:length(CR_ind)
    start_ind = CR_ind(ii);
    end_ind = min(start_ind+25*fps_upsampled, length(NE_seq));
    booster(start_ind:end_ind) = boost_factor;
end
end




