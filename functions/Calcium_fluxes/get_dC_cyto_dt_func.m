%% modelling the partial derivative of c_cyto with respect to time
% return a function handle
function dC_cyto_dt_func = get_dC_cyto_dt_func(theta, other_settings, VGCC_const)
b_ER = theta('b_ER');
b_serca = theta('b_serca');
b_in = theta('b_in');
b_out = theta('b_out');
b_Fmito = theta('b_Fmito');
b_Smito = theta('b_Smito');
K_concentration = theta('K_concentration');
K_serca = theta('K_serca');
K_fmito = theta('K_fmito');
K_smito = theta('K_smito');
b_mPTP = theta('b_mPTP');
ROS_threshold = theta('ROS_threshold');
b_leakage = theta('b_leakage');

FER_func = get_FER_func(theta, K_concentration, other_settings);
Fserca_func = get_Fserca_func(K_serca, other_settings);
F_leakage_func = get_F_leakage_func();
Fmito_func = get_Fmito_func(K_fmito, b_mPTP, ROS_threshold);
Smito_func = get_Smito_func(K_smito);
F_in_func = get_F_in_func(VGCC_const, theta);
F_out_func = get_F_out_func(theta);

dC_cyto_dt_func = @(c_cyto, c_ER, ATP, IP3, PKA_particle, act_cPKC_particle, I_VGCC, AMPA_particle, init_c_ER, c_mito, ROS_mito, ROS_cyto) b_ER * FER_func(c_cyto, c_ER, ATP, IP3, PKA_particle, ROS_cyto) ...
    - b_serca * Fserca_func(c_cyto, c_ER, init_c_ER) ...
    + b_leakage * F_leakage_func(c_cyto, c_ER) ...
    + b_in * F_in_func(I_VGCC, PKA_particle, AMPA_particle) ...
    - b_out * F_out_func(c_cyto, act_cPKC_particle) ...
    + b_Fmito * Fmito_func(c_cyto, c_mito, ROS_mito) ...
    - b_Smito * Smito_func(c_cyto, c_mito);
end