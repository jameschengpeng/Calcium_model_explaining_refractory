%% the function to configure theta
function theta = theta_configuration()
theta = containers.Map;
% param setting, b means coefficient, K means constant in Hill function
theta('b_ER') = 15; % abs value of coeff of FER in dC_cyto_dt
theta('b_serca') = 8; % abs value of coeff of Fsereca in dC_cyto_dt  in literature it is 2.2, but looks like it's too small
theta('b_in') = 1; % abs value of coeff of calcium influx via membrane
theta('b_out') = 0.1; % abs value of coeff of calcium efflux via membrane
theta('b_Fmito') = 2; % abs value of coeff of calcium efflux via mito
theta('b_Smito') = 0.3; % abs value of coeff of calcium influx via mito
theta('K_concentration') = 1; % regulates the effect of concentration diff between c_ER & c_cyto on FER
theta('K_serca') = 0.18; % regulates the serca pump
theta('K_fmito') = 0.2;
theta('K_smito') = 0.2;
theta('K_DA_rew1') = 0.2;
theta('K_DA_rew2') = 3;
theta('b_DA_run') = 0.3;
theta('K_DA_run1') = 0.5;
theta('K_DA_run2') = 2;
theta('b_run') = 0.5;
theta('b_rew') = 1;

theta('K_glu1') = 2;
theta('K_glu2') = 5;
theta('ROS_deletion') = 0.1; % rate of ROS deletion
theta('ROS_delay') = 2;
theta('b_ROS2glu') = 5;
theta('ROS_threshold') = 1; % maximum ROS inside mitochondria
theta('b_ROS_flow') = 1; % coefficient of ROS flow for the change of ROS_mito
theta('b_mPTP') = 1;
theta('b_recep_active') = 0.6;
theta('b_IP3') = 3;
theta('b_IP3_degrade1') = 0.8;
theta('b_IP3_degrade2') = 1.1;
theta('K_IP3_degrade') = 1;
theta('K_IP3_degrade_c_cyto') = 0.3;
% theta('tau_PLC') = 8;
theta('b_cpkc_produce') = 0.5;
theta('b_cpkc_degrade') = 0.9;
theta('b_DAG_degrade') = 10;
theta('b_PLC_produce') = 0.4;
theta('b_PLC_degrade') = 2;
theta('K_IP3_DAG_prod') = 0.2; % Hill affinity value for IP3 and DAG production by PLC
% tau_cpkc_degrade depends on the specific mouse
theta('tau_cpkc_degrade') = 10; % prev 5

theta('c_cyto_peak') = 0.8; % prev 1.2
theta('b_stim2glu_neuron') = 5;
theta('b_glu_neuron2cleft') = 2;
theta('b_glu_cleft2neuron') = 0.05;
theta('tau_cleft2neuron') = 0;
theta('b_glu_cleft2ast') = 0.9;
theta('b_glu_ast2glutamine') = 5;
theta('b_leakage') = 0.005;
% different values of NT delay for different scenarios
theta('NT_delay') = 0;
theta('cPKC_delay') = 1;
theta('NE_thre') = 0;
theta('leaky_factor') = 12; % tau in the leaky integration


theta('b_AC_produce') = 0.8;
theta('b_AC_degrade') = 0.9;
theta('b_cAMP_produce') = 2;
theta('b_cAMP_degrade1') = 0.5;
theta('b_cAMP_degrade2') = 0.5;
theta('b_PKA_produce') = 3.4;
theta('b_PKA_degrade') = 2;
theta('b_PDE_produce') = 0.47;
theta('b_PDE_degrade') = 0.18;
theta('b_GsGPCR_active') = 0.8;

theta('K_AMPA_Glu') = 0.44;

theta('IP3R1_percentage') = 0;
end