%% model the calcium homeostasis in mitochondria
function dC_mito_dt_func = get_dC_mito_dt_func(theta, cyto_mito_ratio, c_mito, c_cyto)
K_fmito = theta('K_fmito');
K_smito = theta('K_smito');
b_Fmito = theta('b_Fmito');
b_Smito = theta('b_Smito');
b_mPTP = theta('b_mPTP');
ROS_threshold = theta('ROS_threshold');

Fmito_func = get_Fmito_func(K_fmito, b_mPTP, c_mito, c_cyto, ROS_threshold);
Smito_func = get_Smito_func(K_smito);

dC_mito_dt_func = @(c_cyto, c_mito, ROS_mito) cyto_mito_ratio * (-b_Fmito * Fmito_func(c_cyto, c_mito, ROS_mito) ...
    + b_Smito * Smito_func(c_cyto, c_mito));
end