%% function to get the time_derivative of m, i.e., the channel activation
function dm_dt_func = get_dm_dt_func()
dm_dt_func = @(m_bar, m, tau_m) (m_bar - m) / tau_m;
end