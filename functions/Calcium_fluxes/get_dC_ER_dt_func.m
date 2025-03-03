%% modelling the partial derivative of c_ER with respect to time
% return a function handle
function dC_ER_dt_func = get_dC_ER_dt_func(theta, cyto_ER_ratio, ROS_cyto, act_cPKC_particle, other_settings)
b_ER = theta('b_ER');
b_serca = theta('b_serca');
K_concentration = theta('K_concentration');
K_serca = theta('K_serca');
b_leakage = theta('b_leakage');

FER_func = get_FER_func(theta, K_concentration, ROS_cyto, act_cPKC_particle, other_settings);
Fserca_func = get_Fserca_func(K_serca, other_settings);
F_leakage_func = get_F_leakage_func();
dC_ER_dt_func = @(c_cyto, c_ER, ATP, IP3, PKA_particle) cyto_ER_ratio * (- b_ER * FER_func(c_cyto, c_ER, ATP, IP3, PKA_particle) ...
    - b_leakage * F_leakage_func(c_cyto, c_ER) ...
    + b_serca * Fserca_func(c_cyto, c_ER));
end