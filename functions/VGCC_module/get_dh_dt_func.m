%% function to get the time_derivative of h, i.e., the channel inactivation
function dh_dt_func = get_dh_dt_func()
dh_dt_func = @(h_bar, h, tau_h) (h_bar - h) / tau_h;
end