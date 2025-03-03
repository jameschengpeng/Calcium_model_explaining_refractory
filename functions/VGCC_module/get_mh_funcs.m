%% get functional forms of m_bar, h_bar, tau_m, tau_h
function [m_T_bar_func, h_T_bar_func, tau_h_Tf_func, tau_h_Ts_func, tau_m_T_func, m_L_bar_func, h_L_func, tau_m_L_func] = get_mh_funcs()
m_T_bar_func = @(V) 1 / (1 + exp(-(V+63.5)/1.5));
h_T_bar_func = @(V) 1 / (1 + exp((V+76.2)/3));
tau_h_Tf_func = @(V) 50 * exp(-((V+72)/10)^2) + 10;
tau_h_Ts_func = @(V) 400 * exp(-((V+100)/10)^2) + 400;
tau_m_T_func = @(V) 65 * exp(-((V+68)/6)^2) + 12;
m_L_bar_func = @(V) 1 / (1 + exp(-(V+50)/3));
h_L_func = @(c_cyto) 0.00045 / (0.00045 + c_cyto);
tau_m_L_func = @(V) 18 * exp(-((V+45)/20)^2) + 1.5;
end