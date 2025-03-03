%% the function that executes an entire episode
% the numerical solution is computed by the 4th order Runge-Kutta method
function [chemical_table, flux_table] = run(theta, other_settings, run_seq, NE_spont, running_threshold, fps_org, fps_upsampled, ...
    stim_onset, cyto_ER_ratio, cyto_mito_ratio, init_IP3, DA_spont, reward_onset, ...
    init_c_cyto, init_c_ER, init_c_mito, init_ROS_mito, init_ROS_cyto, IP3_spont, smoothing_factor, TT)
%% interpretation of the parameters
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
K_DA_rew1 = theta('K_DA_rew1');
K_DA_rew2 = theta('K_DA_rew2');
b_DA_run = theta('b_DA_run');
K_DA_run1 = theta('K_DA_run1');
K_DA_run2 = theta('K_DA_run2');
b_run = theta('b_run');
b_rew = theta('b_rew');
b_IP3 = theta('b_IP3');
K_glu1 = theta('K_glu1');
K_glu2 = theta('K_glu2');
ROS_deletion = theta('ROS_deletion'); % rate of ROS deletion
ROS_delay = theta('ROS_delay');
b_ROS2glu = theta('b_ROS2glu');
ROS_threshold = theta('ROS_threshold'); % maximum ROS inside mitochondria
b_ROS_flow = theta('b_ROS_flow');
b_mPTP = theta('b_mPTP');
b_recep_active = theta('b_recep_active');
b_IP3_degrade1 = theta('b_IP3_degrade1');
b_IP3_degrade2 = theta('b_IP3_degrade2');
K_IP3_degrade = theta('K_IP3_degrade');
K_IP3_degrade_c_cyto = theta('K_IP3_degrade_c_cyto');
% tau_PLC = theta('tau_PLC');
b_cpkc_produce = theta('b_cpkc_produce');
b_cpkc_degrade = theta('b_cpkc_degrade');
b_DAG_degrade = theta('b_DAG_degrade');
b_PLC_produce = theta('b_PLC_produce');
b_PLC_degrade = theta('b_PLC_degrade');
K_IP3_DAG_prod = theta('K_IP3_DAG_prod');
tau_cpkc_degrade = theta('tau_cpkc_degrade');
c_cyto_peak = theta('c_cyto_peak');
b_stim2glu_neuron = theta('b_stim2glu_neuron');
b_glu_neuron2cleft = theta('b_glu_neuron2cleft');
b_glu_cleft2neuron = theta('b_glu_cleft2neuron');
tau_cleft2neuron = theta('tau_cleft2neuron');
b_glu_cleft2ast = theta('b_glu_cleft2ast');
b_glu_ast2glutamine = theta('b_glu_ast2glutamine');
b_leakage = theta('b_leakage');
b_AC_produce = theta('b_AC_produce');
b_AC_degrade = theta('b_AC_degrade');
b_cAMP_produce = theta('b_cAMP_produce');
b_cAMP_degrade1 = theta('b_cAMP_degrade1');
b_cAMP_degrade2 = theta('b_cAMP_degrade2');
b_PKA_produce = theta('b_PKA_produce');
b_PKA_degrade = theta('b_PKA_degrade');
b_PDE_produce = theta('b_PDE_produce');
b_PDE_degrade = theta('b_PDE_degrade');
NT_delay = theta('NT_delay');
cPKC_delay = theta('cPKC_delay');
NE_thre = theta('NE_thre');
leaky_factor = theta('leaky_factor');
b_GsGPCR_active = theta('b_GsGPCR_active');
K_AMPA_Glu = theta('K_AMPA_Glu');


VGCC_const = get_VGCC_const();
%%
dt = 1/fps_upsampled;
run_func = seq2run_func(run_seq, fps_org);
% [NE_seq, NE_func] = run2NE_seq_func(run_seq, NE_spont, running_threshold, fps_upsampled, smoothing_factor);
[NE_seq, NE_func] = run2NE_seq_func(run_seq, NE_spont, fps_org, fps_upsampled, TT, theta);
[glu_stim_seq, glu_stim_func] = run2glu_stim_seq_func(run_seq, stim_onset, fps_org, fps_upsampled, K_glu1, K_glu2, other_settings, theta);
[DA_seq, DA_func] = run2DA_seq_func(run_seq, DA_spont, reward_onset, fps_org, fps_upsampled, running_threshold, smoothing_factor, ...
    K_DA_rew1, K_DA_rew2, b_DA_run, b_run, b_rew);
termination_time = (length(run_seq)-1)/fps_org;
n_steps = ceil(fps_upsampled * termination_time);
discrete_time = linspace(0, termination_time, n_steps);
run_seq = arrayfun(run_func, discrete_time); % upsampled run_seq

IP3 = init_IP3;
c_cyto = init_c_cyto;
c_ER = init_c_ER;
c_mito = init_c_mito;
ROS_mito = init_ROS_mito;
ROS_cyto = init_ROS_cyto;
R_particle = [ones(100,1) zeros(100,1)];
PLC_particle = [ones(500,1) zeros(500,1)]; % initially all PLC particles are inactive
act_cPKC_particle = [ones(500,1) zeros(500,1)]; % initially all cPKC particles are inactive;
DAG = 0;

AMPA_particle = [ones(100,1) zeros(100,1) zeros(100,1)]; % AMPA receptor has three states: desensitized close open

glu_neuron = 0;
glu_cleft = 0;
glu_ast = 0;
glutamine_ast = 0;

AC_particle = [ones(100,1) zeros(100,1)];
PKA_particle = [ones(100,1) zeros(100,1)];
GsGPCR_particle = [ones(100,1) zeros(100,1)];
PDE_particle = [ones(100,1) zeros(100,1)];
cAMP = 0;

% VGCC components initialization
V = -70;
E_Ca_func = get_E_Ca_func(VGCC_const); E_Ca = E_Ca_func(c_cyto);
m_T = 0; % m_T = 1 / (1 + exp(-(V+63.5)/1.5));
h_Tf = 0; % h_Tf = 1 / (1 + exp((V+76.2)/3));
h_Ts = 0; % h_Ts = 1 / (1 + exp((V+76.2)/3));
m_L = 0; % m_L = 1 / (1 + exp(-(V+50)/3));
h_L = 0; % h_L = 0.00045 / (0.00045 + c_cyto);

I_T_type_func = get_I_T_type_func(VGCC_const); I_T_type = I_T_type_func(V, E_Ca, m_T, h_Tf, h_Ts);
I_L_type_func = get_I_L_type_func(VGCC_const); I_L_type = I_L_type_func(V, E_Ca, m_L, h_L);
I_VGCC = I_L_type + I_T_type;
%% chemical substance
c_cyto_record = zeros(size(run_seq)); c_cyto_record(1) = c_cyto;
c_ER_record = zeros(size(run_seq)); c_ER_record(1) = c_ER;
c_mito_record = zeros(size(run_seq)); c_mito_record(1) = c_mito;
IP3_record = zeros(size(run_seq)); IP3_record(1) = IP3;
ATP_func = get_ATP_func(NE_func, glu_ast);
ATP = ATP_func(0);
ATP_record = zeros(size(run_seq)); ATP_record(1) = ATP;
ROS_mito_record = zeros(size(run_seq)); ROS_mito_record(1) = ROS_mito;
ROS_cyto_record = zeros(size(run_seq)); ROS_cyto_record(1) = ROS_cyto;
R_record = zeros(size(run_seq)); R_record(1) = mean(R_particle(:,2));
PLC_record = zeros(size(run_seq)); PLC_record(1) = mean(PLC_particle(:,2));
act_cPKC_record = zeros(size(run_seq)); act_cPKC_record(1) = mean(act_cPKC_particle(:,2));
DAG_record = zeros(size(run_seq)); DAG_record(1) = DAG;

AMPA_record = zeros(size(run_seq)); AMPA_record(1) = mean(AMPA_particle(:,3));

glu_neuron_record = zeros(size(run_seq)); glu_neuron_record(1) = glu_neuron;
glu_cleft_record = zeros(size(run_seq)); glu_cleft_record(1) = glu_cleft;
glu_ast_record = zeros(size(run_seq)); glu_ast_record(1) = glu_ast;
glutamine_ast_record = zeros(size(run_seq)); glutamine_ast_record(1) = glutamine_ast;

AC_record = zeros(size(run_seq)); AC_record(1) = mean(AC_particle(:,2));
PKA_record = zeros(size(run_seq)); PKA_record(1) = mean(PKA_particle(:,2));
GsGPCR_record = zeros(size(run_seq)); GsGPCR_record(1) = mean(GsGPCR_particle(:,2));
PDE_record = zeros(size(run_seq)); PDE_record(1) = mean(PDE_particle(:,2));
cAMP_record = zeros(size(run_seq)); cAMP_record(1) = cAMP;

V_record = zeros(size(run_seq)); V_record(1) = V;
I_VGCC_record = zeros(size(run_seq)); I_VGCC_record(1) = I_VGCC;
%% flux
FER_record = zeros(size(run_seq));
Fserca_record = zeros(size(run_seq));
F_leakage_record = zeros(size(run_seq));
Fmito_record = zeros(size(run_seq));
Smito_record = zeros(size(run_seq));
F_in_record = zeros(size(run_seq));
F_out_record = zeros(size(run_seq));
c_cyto_eff_record = zeros(size(run_seq));
ATP_eff_record = zeros(size(run_seq));
IP3_eff_record = zeros(size(run_seq));

for ii = 2:length(NE_seq)
    NE = NE_seq(ii); glu_stim = glu_stim_seq(ii); DA = DA_seq(ii);
    t0 = (ii-2) * dt; % the existing time point
    %% update on glu_neuron: d_glu_neuron_dt_func = @(glu_neuron, glu_cleft_record, t, fps_upsampled)
    d_glu_neuron_dt_func = get_d_glu_neuron_dt_func(theta, glu_stim_func);
    %% update on glu_cleft: d_glu_cleft_dt_func = @(glu_cleft_record, glu_neuron, t, fps_upsampled)
    d_glu_cleft_dt_func = get_d_glu_cleft_dt_func(theta);
    %% update on glu_ast: d_glu_ast_dt_func = @(glu_ast, glu_cleft)
    d_glu_ast_dt_func = get_d_glu_ast_dt_func(theta);
    %% update on glutamine_ast: d_glutamine_ast_dt_func = @(glu_ast, glutamine_ast)
    d_glutamine_ast_dt_func = get_d_glutamine_ast_dt_func(theta);

    k1_glu_neuron = dt * d_glu_neuron_dt_func(glu_neuron, glu_cleft_record, t0, fps_upsampled);
    k1_glu_cleft = dt * d_glu_cleft_dt_func(glu_cleft_record, glu_cleft, glu_neuron, t0, fps_upsampled);
    k1_glu_ast = dt * d_glu_ast_dt_func(glu_ast, glu_cleft);
    k1_glutamine_ast = dt * d_glutamine_ast_dt_func(glu_ast, glutamine_ast);
    
    k2_glu_neuron = dt * d_glu_neuron_dt_func(glu_neuron+k1_glu_neuron/2, glu_cleft_record, t0+dt/2, fps_upsampled);
    k2_glu_cleft = dt * d_glu_cleft_dt_func(glu_cleft_record, glu_cleft+k1_glu_cleft/2, glu_neuron+k1_glu_neuron/2, t0+dt/2, fps_upsampled);
    k2_glu_ast = dt * d_glu_ast_dt_func(glu_ast+k1_glu_ast/2, glu_cleft+k1_glu_cleft/2);
    k2_glutamine_ast = dt * d_glutamine_ast_dt_func(glu_ast+k1_glu_ast/2, glutamine_ast+k1_glutamine_ast/2);

    k3_glu_neuron = dt * d_glu_neuron_dt_func(glu_neuron+k2_glu_neuron/2, glu_cleft_record, t0+dt/2, fps_upsampled);
    k3_glu_cleft = dt * d_glu_cleft_dt_func(glu_cleft_record, glu_cleft+k2_glu_cleft/2, glu_neuron+k2_glu_neuron/2, t0+dt/2, fps_upsampled);
    k3_glu_ast = dt * d_glu_ast_dt_func(glu_ast+k2_glu_ast/2, glu_cleft+k2_glu_cleft/2);
    k3_glutamine_ast = dt * d_glutamine_ast_dt_func(glu_ast+k2_glu_ast/2, glutamine_ast+k2_glutamine_ast/2);

    k4_glu_neuron = dt * d_glu_neuron_dt_func(glu_neuron+k3_glu_neuron, glu_cleft_record, t0+dt, fps_upsampled);
    k4_glu_cleft = dt * d_glu_cleft_dt_func(glu_cleft_record, glu_cleft+k3_glu_cleft, glu_neuron+k3_glu_neuron, t0+dt, fps_upsampled);
    k4_glu_ast = dt * d_glu_ast_dt_func(glu_ast+k3_glu_ast, glu_cleft+k3_glu_cleft);
    k4_glutamine_ast = dt * d_glutamine_ast_dt_func(glu_ast+k3_glu_ast, glutamine_ast+k3_glutamine_ast);

    glu_neuron_new = max(0.01, glu_neuron + k1_glu_neuron/6 + k2_glu_neuron/3 + k3_glu_neuron/3 + k4_glu_neuron/6);
    glu_cleft_new = max(0.01, glu_cleft + k1_glu_cleft/6 + k2_glu_cleft/3 + k3_glu_cleft/3 + k4_glu_cleft/6);
    glu_ast_new = max(0.01, glu_ast + k1_glu_ast/6 + k2_glu_ast/3 + k3_glu_ast/3 + k4_glu_ast/6);
    glutamine_ast_new = max(0.01, glutamine_ast + k1_glutamine_ast/6 + k2_glutamine_ast/3 + k3_glutamine_ast/3 + k4_glutamine_ast/6);

    %% update on R_particle (activated % of Gq-GPCR), PLC_particle, IP3; dR_particle_dt_func = @(act_cPKC_particle, R_particle, t); 
    %% dPLC_particle_dt_func = @(act_cPKC_particle, PLC_particle, R_particle); dIP3_dt_func = @(c_cyto, IP3, PLC_particle)
    %% update on activated cPKC (proportion), DAG. d_act_cPKC_particle_dt_func = @(c_cyto, c_cyto_record, DAG, DAG_record, act_cPKC_particle, fps_upsampled, ii); dDAG_dt_func = @(c_cyto, PLC_particle, DAG)
    %% update on AMPA; d_AMPA_particle_dt_func = @(AMPA_particle, glu_cleft)
    dR_particle_dt_func = get_dR_particle_dt_func(theta, other_settings, NE_func, NE_spont, DA_func, DA_spont, glu_cleft, t0);
    dPLC_particle_dt_func = get_dPLC_particle_dt_func(theta);
    dIP3_dt_func = get_dIP3_dt_func(theta);
    d_act_cPKC_particle_dt_func = get_d_act_cPKC_particle_dt_func(theta);
    dDAG_dt_func = get_dDAG_dt_func(theta);
    d_AMPA_particle_dt_func = get_d_AMPA_particle_dt_func(theta);

    k1_R = dt * dR_particle_dt_func(act_cPKC_particle, R_particle, t0);
    k1_PLC = dt * dPLC_particle_dt_func(act_cPKC_particle, PLC_particle, R_particle);
    k1_IP3 = dt * dIP3_dt_func(c_cyto, IP3, PLC_particle);
    k1_cPKC = dt * d_act_cPKC_particle_dt_func(c_cyto, c_cyto_record, DAG, DAG_record, act_cPKC_particle, fps_upsampled, ii);
    k1_DAG = dt * dDAG_dt_func(c_cyto, PLC_particle, DAG);
    k1_AMPA = dt * d_AMPA_particle_dt_func(AMPA_particle, glu_cleft);

    k2_PLC = dt * dPLC_particle_dt_func(act_cPKC_particle+k1_cPKC/2, PLC_particle+k1_PLC/2, R_particle+k1_R/2);
    k2_R = dt * dR_particle_dt_func(act_cPKC_particle+k1_cPKC/2, R_particle+k1_R/2, t0+dt/2);
    k2_IP3 = dt * dIP3_dt_func(c_cyto, IP3 + 1/2*k1_IP3, PLC_particle + 1/2*k1_PLC);
    k2_cPKC = dt * d_act_cPKC_particle_dt_func(c_cyto, c_cyto_record, DAG+k1_DAG/2, DAG_record, act_cPKC_particle+k1_cPKC/2, fps_upsampled, ii);
    k2_DAG = dt * dDAG_dt_func(c_cyto, PLC_particle+k1_PLC/2, DAG+k1_DAG/2);
    k2_AMPA = dt * d_AMPA_particle_dt_func(AMPA_particle+k1_AMPA/2, glu_cleft+k1_glu_cleft/2);

    k3_R = dt * dR_particle_dt_func(act_cPKC_particle+k2_cPKC/2, R_particle+k2_R/2, t0+dt/2);
    k3_PLC = dt * dPLC_particle_dt_func(act_cPKC_particle+k2_cPKC/2, PLC_particle+k2_PLC/2, R_particle+k2_R/2);
    k3_IP3 = dt * dIP3_dt_func(c_cyto, IP3 + 1/2*k2_IP3, PLC_particle + 1/2*k2_PLC);
    k3_cPKC = dt * d_act_cPKC_particle_dt_func(c_cyto, c_cyto_record, DAG+k2_DAG/2, DAG_record, act_cPKC_particle+k2_cPKC/2, fps_upsampled, ii);
    k3_DAG = dt * dDAG_dt_func(c_cyto, PLC_particle+k2_PLC/2, DAG+k2_DAG/2);
    k3_AMPA = dt * d_AMPA_particle_dt_func(AMPA_particle+k2_AMPA/2, glu_cleft+k2_glu_cleft/2);

    k4_R = dt * dR_particle_dt_func(act_cPKC_particle+k3_cPKC, R_particle+k3_R, t0+dt);
    k4_PLC = dt * dPLC_particle_dt_func(act_cPKC_particle+k3_cPKC, PLC_particle+k3_PLC, R_particle+k3_R);
    k4_IP3 = dt * dIP3_dt_func(c_cyto, IP3 + k3_IP3, PLC_particle+k3_PLC);
    k4_cPKC = dt * d_act_cPKC_particle_dt_func(c_cyto, c_cyto_record, DAG+k3_DAG, DAG_record, act_cPKC_particle+k3_cPKC, fps_upsampled, ii);
    k4_DAG = dt * dDAG_dt_func(c_cyto, PLC_particle+k3_PLC, DAG+k3_DAG);
    k4_AMPA = dt * d_AMPA_particle_dt_func(AMPA_particle+k3_AMPA, glu_cleft+k3_glu_cleft);
    
    R_new = R_particle + k1_R/6 + k2_R/3 + k3_R/3 + k4_R/6;
    PLC_new = PLC_particle + k1_PLC/6 + k2_PLC/3 + k3_PLC/3 + k4_PLC/6;
    IP3_new = max(IP3_spont, IP3 + k1_IP3/6 + k2_IP3/3 + k3_IP3/3 + k4_IP3/6);
    act_cPKC_new = act_cPKC_particle + k1_cPKC/6 + k2_cPKC/3 + k3_cPKC/3 + k4_cPKC/6;
    DAG_new = max(0, DAG + k1_DAG/6 + k2_DAG/3 + k3_DAG/3 + k4_DAG/6);
    AMPA_new = AMPA_particle + k1_AMPA/6 + k2_AMPA/3 + k3_AMPA/3 + k4_AMPA/6;

    %% update on AC, cAMP, PKA, PDE, GsGPCR
    % dAC_particle_dt_func = @(PKA_particle, AC_particle, R_particle); dcAMP_dt_func = @(cAMP, AC_particle, PDE_particle)
    % dPKA_particle_dt_func = @(PKA_particle, cAMP);
    % dPDE_particle_dt_func = @(PDE_particle, c_cyto); dGsGPCR_particle_dt_func = @(PKA_particle, GsGPCR_particle, t)
    dGsGPCR_particle_dt_func = get_dGsGPCR_particle_dt_func(theta, NE_func, NE_spont, DA_func, DA_spont, t0);
    dAC_particle_dt_func = get_dAC_particle_dt_func(theta);
    dcAMP_dt_func = get_dcAMP_dt_func(theta);
    dPKA_particle_dt_func = get_dPKA_particle_dt_func(theta);
    dPDE_particle_dt_func = get_dPDE_particle_dt_func(theta);

    k1_GsGPCR = dt * dGsGPCR_particle_dt_func(PKA_particle, GsGPCR_particle, t0);
    k1_AC = dt * dAC_particle_dt_func(PKA_particle, AC_particle, GsGPCR_particle);
    k1_cAMP = dt * dcAMP_dt_func(cAMP, AC_particle, PDE_particle);
    k1_PKA = dt * dPKA_particle_dt_func(PKA_particle, cAMP);
    k1_PDE = dt * dPDE_particle_dt_func(PDE_particle, c_cyto);

    k2_GsGPCR = dt * dGsGPCR_particle_dt_func(PKA_particle+k1_PKA/2, GsGPCR_particle+k1_GsGPCR/2, t0+dt/2);
    k2_AC = dt * dAC_particle_dt_func(PKA_particle+k1_PKA/2, AC_particle+k1_AC/2, GsGPCR_particle+k1_GsGPCR/2);
    k2_cAMP = dt * dcAMP_dt_func(cAMP+k1_cAMP/2, AC_particle+k1_AC/2, PDE_particle+k1_PDE/2);
    k2_PKA = dt * dPKA_particle_dt_func(PKA_particle+k1_PKA/2, cAMP+k1_cAMP/2);
    k2_PDE = dt * dPDE_particle_dt_func(PDE_particle+k1_PDE/2, c_cyto);

    k3_GsGPCR = dt * dGsGPCR_particle_dt_func(PKA_particle+k2_PKA/2, GsGPCR_particle+k2_GsGPCR/2, t0+dt/2);
    k3_AC = dt * dAC_particle_dt_func(PKA_particle+k2_PKA/2, AC_particle+k2_AC/2, GsGPCR_particle+k2_GsGPCR/2);
    k3_cAMP = dt * dcAMP_dt_func(cAMP+k2_cAMP/2, AC_particle+k2_AC/2, PDE_particle+k2_PDE/2);
    k3_PKA = dt * dPKA_particle_dt_func(PKA_particle+k2_PKA/2, cAMP+k2_cAMP/2);
    k3_PDE = dt * dPDE_particle_dt_func(PDE_particle+k2_PDE/2, c_cyto);

    k4_GsGPCR = dt * dGsGPCR_particle_dt_func(PKA_particle+k3_PKA, GsGPCR_particle+k3_GsGPCR, t0+dt);
    k4_AC = dt * dAC_particle_dt_func(PKA_particle+k3_PKA, AC_particle+k3_AC, GsGPCR_particle+k3_GsGPCR);
    k4_cAMP = dt * dcAMP_dt_func(cAMP+k3_cAMP, AC_particle+k3_AC, PDE_particle+k3_PDE);
    k4_PKA = dt * dPKA_particle_dt_func(PKA_particle+k3_PKA, cAMP+k3_cAMP);
    k4_PDE = dt * dPDE_particle_dt_func(PDE_particle+k3_PDE, c_cyto);


    GsGPCR_new = GsGPCR_particle + k1_GsGPCR/6 + k2_GsGPCR/3 + k3_GsGPCR/3 + k4_GsGPCR/6;
    AC_new = AC_particle + k1_AC/6 + k2_AC/3 + k3_AC/3 + k4_AC/6;
    cAMP_new = cAMP + k1_cAMP/6 + k2_cAMP/3+ k3_cAMP/3 + k4_cAMP/6;
    PKA_new = PKA_particle + k1_PKA/6 + k2_PKA/3 + k3_PKA/3 + k4_PKA/6;
    PDE_new = PDE_particle + k1_PDE/6 + k2_PDE/3 + k3_PDE/3 + k4_PDE/6;

    %% update on ATP, ATP_func = @(t)
    ATP_func = get_ATP_func(NE_func, glu_ast);
    ATP_new = ATP_func(t0 + dt);

    %% update on the VGCC module
    % update the tau and steady state values
    [m_T_bar_func, h_T_bar_func, tau_h_Tf_func, tau_h_Ts_func, tau_m_T_func, m_L_bar_func, h_L_func, tau_m_L_func] = get_mh_funcs();
    dm_dt_func = get_dm_dt_func(); % dm_dt_func = @(m_bar, m, tau_m)
    dh_dt_func = get_dh_dt_func(); % dh_dt_func = @(h_bar, h, tau_h)
    E_Ca_func = get_E_Ca_func(VGCC_const); % E_Ca_func = @(c_cyto)
    dV_dt_func = get_dV_dt_func(VGCC_const); % dV_dt_func = @(I_VGCC)
    I_VGCC_func = get_I_VGCC_func(); % I_VGCC_func = @(I_T_type, I_L_type)
    I_T_type_func = get_I_T_type_func(VGCC_const); % I_T_type_func = @(V, E_Ca, m_T, h_Tf, h_Ts)
    I_L_type_func = get_I_L_type_func(VGCC_const); % I_L_type_func = @(V, E_Ca, m_L, h_L)


    I_T_type1 = I_T_type_func(V, E_Ca, m_T, h_Tf, h_Ts);
    I_L_type1 = I_L_type_func(V, E_Ca, m_L, h_L);
    k1_V = dt * dV_dt_func(I_VGCC_func(I_T_type1, I_L_type1));
    k1_m_T = dt * dm_dt_func(m_T_bar_func(V), m_T, tau_m_T_func(V));
    k1_h_Tf = dt * dh_dt_func(h_T_bar_func(V), h_Tf, tau_h_Tf_func(V));
    k1_h_Ts = dt * dh_dt_func(h_T_bar_func(V), h_Ts, tau_h_Ts_func(V));
    k1_m_L = dt * dm_dt_func(m_L_bar_func(V), m_L, tau_m_L_func(V));

    I_T_type2 = I_T_type_func(V + k1_V/2, E_Ca, m_T + k1_m_T/2, h_Tf + k1_h_Tf/2, h_Ts + k1_h_Ts/2);
    I_L_type2 = I_L_type_func(V + k1_V/2, E_Ca, m_L + k1_m_L/2, h_L);
    k2_V = dt * dV_dt_func(I_VGCC_func(I_T_type2, I_L_type2));
    k2_m_T = dt * dm_dt_func(m_T_bar_func(V + k1_V/2), m_T + k1_m_T/2, tau_m_T_func(V + k1_V/2));
    k2_h_Tf = dt * dh_dt_func(h_T_bar_func(V + k1_V/2), h_Tf + k1_h_Tf/2, tau_h_Tf_func(V + k1_V/2));
    k2_h_Ts = dt * dh_dt_func(h_T_bar_func(V + k1_V/2), h_Ts + k1_h_Ts/2, tau_h_Ts_func(V + k1_V/2));
    k2_m_L = dt * dm_dt_func(m_L_bar_func(V + k1_V/2), m_L + k1_m_L/2, tau_m_L_func(V + k1_V/2));

    I_T_type3 = I_T_type_func(V + k2_V/2, E_Ca, m_T + k2_m_T/2, h_Tf + k2_h_Tf/2, h_Ts + k2_h_Ts/2);
    I_L_type3 = I_L_type_func(V + k2_V/2, E_Ca, m_L + k2_m_L/2, h_L);
    k3_V = dt * dV_dt_func(I_VGCC_func(I_T_type3, I_L_type3));
    k3_m_T = dt * dm_dt_func(m_T_bar_func(V + k2_V/2), m_T + k2_m_T/2, tau_m_T_func(V + k2_V/2));
    k3_h_Tf = dt * dh_dt_func(h_T_bar_func(V + k2_V/2), h_Tf + k2_h_Tf/2, tau_h_Tf_func(V + k2_V/2));
    k3_h_Ts = dt * dh_dt_func(h_T_bar_func(V + k2_V/2), h_Ts + k2_h_Ts/2, tau_h_Ts_func(V + k2_V/2));
    k3_m_L = dt * dm_dt_func(m_L_bar_func(V + k2_V/2), m_L + k2_m_L/2, tau_m_L_func(V + k2_V/2));

    I_T_type4 = I_T_type_func(V + k3_V, E_Ca, m_T + k3_m_T, h_Tf + k3_h_Tf, h_Ts + k3_h_Ts);
    I_L_type4 = I_L_type_func(V + k3_V, E_Ca, m_L + k3_m_L, h_L);
    k4_V = dt * dV_dt_func(I_VGCC_func(I_T_type4, I_L_type4));
    k4_m_T = dt * dm_dt_func(m_T_bar_func(V + k3_V), m_T + k3_m_T, tau_m_T_func(V + k3_V));
    k4_h_Tf = dt * dh_dt_func(h_T_bar_func(V + k3_V), h_Tf + k3_h_Tf, tau_h_Tf_func(V + k3_V));
    k4_h_Ts = dt * dh_dt_func(h_T_bar_func(V + k3_V), h_Ts + k3_h_Ts, tau_h_Ts_func(V + k3_V));
    k4_m_L = dt * dm_dt_func(m_L_bar_func(V + k3_V), m_L + k3_m_L, tau_m_L_func(V + k3_V));

    V_new = min(-65, V + k1_V/6 + k2_V/3 + k3_V/3 + k4_V/6);
    % V_new = V + k1_V/6 + k2_V/3 + k3_V/3 + k4_V/6;
    m_T_new = m_T + k1_m_T/6 + k2_m_T/3 + k3_m_T/3 + k4_m_T/6;
    h_Tf_new = h_Tf + k1_h_Tf/6 + k2_h_Tf/3 + k3_h_Tf/3 + k4_h_Tf/6;
    h_Ts_new = h_Ts + k1_h_Ts/6 + k2_h_Ts/3 + k3_h_Ts/3 + k4_h_Ts/6;
    m_L_new = m_L + k1_m_L/6 + k2_m_L/3 + k3_m_L/3 + k4_m_L/6;

    I_VGCC = I_VGCC_func(I_T_type1, I_L_type1);
    % update of E_Ca, h_L, I_T_type, I_L_type, I_VGCC are after c_cyto update

    %% update on FER, FER_func = @(c_cyto, c_ER, ATP, IP3, PKA_particle)
    [FER_func, c_cyto_effect, ATP_effect, IP3_effect] = get_FER_func(theta, K_concentration, ROS_cyto, act_cPKC_particle, other_settings);
    c_cyto_effect = c_cyto_effect(c_cyto);
    ATP_effect = ATP_effect(ATP, IP3);
    FER = FER_func(c_cyto, c_ER, ATP, IP3, PKA_particle);

    %% update on Fserca, Fserca_func = @(c_cyto)
    Fserca_func = get_Fserca_func(K_serca, other_settings);
    Fserca = Fserca_func(c_cyto, c_ER);

    %% update on F_leakage, F_leakage_func = @(c_cyto, c_ER)
    F_leakage_func = get_F_leakage_func();
    F_leakage = F_leakage_func(c_cyto, c_ER);

    %% update on F_in and F_out  F_in_func = @(IP3)    F_out_func = @(c_cyto)
    F_in_func = get_F_in_func(VGCC_const);
    F_out_func = get_F_out_func();
    F_in = F_in_func(I_VGCC, PKA_particle, AMPA_particle);
    F_out = F_out_func(c_cyto, act_cPKC_particle);

    %% update on Fmito, Fmito_func = @(c_cyto, c_mito, ROS_mito)
    Fmito_func = get_Fmito_func(K_fmito, b_mPTP, c_mito, c_cyto, ROS_threshold);
    Fmito = Fmito_func(c_cyto, c_mito, ROS_mito);

    %% update on Smito, Smito_func = @(c_cyto)
    Smito_func = get_Smito_func(K_smito);
    Smito = Smito_func(c_cyto, c_mito);

    %% update on ROS_mito, ROS_cyto; dROS_cyto_dt_func = @(ROS_mito, ROS_cyto); dROS_mito_dt_func = @(ROS_mito, ROS_cyto, c_mito, glu_ast_record, t)
    dROS_mito_dt_func = get_dROS_mito_dt_func(theta, fps_upsampled);
    dROS_cyto_dt_func = get_dROS_cyto_dt_func(theta, cyto_mito_ratio);

    %% update on c_cyto, c_ER, and c_mito, dC_cyto_dt_func = @(c_cyto, c_ER, ATP, IP3); dC_ER_dt_func = @(c_cyto, c_ER, ATP, IP3); dC_mito_dt_func = @(c_cyto, c_mito, ROS_mito)
    dC_cyto_dt_func = get_dC_cyto_dt_func(theta, c_mito, c_cyto, ROS_mito, ROS_cyto, act_cPKC_particle, other_settings, VGCC_const);
    dC_ER_dt_func = get_dC_ER_dt_func(theta, cyto_ER_ratio, ROS_cyto, act_cPKC_particle, other_settings);
    dC_mito_dt_func = get_dC_mito_dt_func(theta, cyto_mito_ratio, c_mito, c_cyto);

    k1_ROS_mito = dt * dROS_mito_dt_func(ROS_mito, ROS_cyto, c_mito, glu_ast_record, t0);
    k1_ROS_cyto = dt * dROS_cyto_dt_func(ROS_mito, ROS_cyto);
    k1_c_cyto = dt * dC_cyto_dt_func(c_cyto, c_ER, ATP_func(t0), IP3, PKA_particle, act_cPKC_particle, I_VGCC, AMPA_particle);
    k1_c_ER = dt * dC_ER_dt_func(c_cyto, c_ER, ATP_func(t0), IP3, PKA_particle);
    k1_c_mito = dt * dC_mito_dt_func(c_cyto, c_mito, ROS_mito);

    k2_ROS_mito = dt * dROS_mito_dt_func(ROS_mito + k1_ROS_mito/2, ROS_cyto + k1_ROS_cyto/2, c_mito + k1_c_mito/2, glu_ast_record, t0+dt/2);
    k2_ROS_cyto = dt * dROS_cyto_dt_func(ROS_mito + k1_ROS_mito/2, ROS_cyto + k1_ROS_cyto/2);
    k2_c_cyto = dt * dC_cyto_dt_func(max(0,c_cyto + k1_c_cyto/2), max(0,c_ER + k1_c_ER/2), ATP_func(t0 + dt/2), IP3+k1_IP3/2, PKA_particle+k1_PKA/2, act_cPKC_particle+k1_cPKC/2, I_VGCC, AMPA_particle+k1_AMPA/2);
    k2_c_ER = dt * dC_ER_dt_func(max(0,c_cyto + k1_c_cyto/2), max(0,c_ER + k1_c_ER/2), ATP_func(t0 + dt/2), IP3+k1_IP3/2, PKA_particle+k1_PKA/2);
    k2_c_mito = dt * dC_mito_dt_func(max(0,c_cyto+k1_c_cyto/2), max(0,c_mito+k1_c_mito/2), ROS_mito+k1_ROS_mito/2);

    k3_ROS_mito = dt * dROS_mito_dt_func(ROS_mito + k2_ROS_mito/2, ROS_cyto + k2_ROS_cyto/2, c_mito + k2_c_mito/2, glu_ast_record, t0+dt/2);
    k3_ROS_cyto = dt * dROS_cyto_dt_func(ROS_mito + k2_ROS_mito/2, ROS_cyto + k2_ROS_cyto/2);
    k3_c_cyto = dt * dC_cyto_dt_func(max(0,c_cyto + k2_c_cyto/2), max(0,c_ER + k2_c_ER/2), ATP_func(t0 + dt/2), IP3+k2_IP3/2, PKA_particle+k2_PKA/2, act_cPKC_particle+k2_cPKC/2, I_VGCC, AMPA_particle+k2_AMPA/2);
    k3_c_ER = dt * dC_ER_dt_func(max(0,c_cyto + k2_c_cyto/2), max(0,c_ER + k2_c_ER/2), ATP_func(t0 + dt/2), IP3+k2_IP3/2, PKA_particle+k2_PKA/2);
    k3_c_mito = dt * dC_mito_dt_func(max(0,c_cyto+k2_c_cyto/2), max(0,c_mito+k2_c_mito/2), ROS_mito+k2_ROS_mito/2);

    k4_ROS_mito = dt * dROS_mito_dt_func(ROS_mito + k3_ROS_mito, ROS_cyto + k3_ROS_cyto, c_mito + k3_c_mito, glu_ast_record, t0+dt);
    k4_ROS_cyto = dt * dROS_cyto_dt_func(ROS_mito + k3_ROS_mito, ROS_cyto + k3_ROS_cyto);
    k4_c_cyto = dt * dC_cyto_dt_func(max(0,c_cyto + k3_c_cyto), max(0,c_ER + k3_c_ER), ATP_func(t0 + dt), IP3+k3_IP3, PKA_particle+k3_PKA, act_cPKC_particle+k3_cPKC, I_VGCC, AMPA_particle+k3_AMPA);
    k4_c_ER = dt * dC_ER_dt_func(max(0,c_cyto + k3_c_cyto), max(0,c_ER + k3_c_ER), ATP_func(t0 + dt), IP3+k3_IP3, PKA_particle+k3_PKA);
    k4_c_mito = dt * dC_mito_dt_func(max(0,c_cyto+k3_c_cyto), max(0,c_mito+k3_c_mito), ROS_mito+k3_ROS_mito);


    ROS_mito_new = max(0, ROS_mito + k1_ROS_mito/6 + k2_ROS_mito/3 + k3_ROS_mito/3 + k4_ROS_mito/6);
    ROS_cyto_new = max(0, ROS_cyto + k1_ROS_cyto/6 + k2_ROS_cyto/3 + k3_ROS_cyto/3 + k4_ROS_cyto/6);

    cyto_change = sum([k1_c_cyto/6, k2_c_cyto/3, k3_c_cyto/3, k4_c_cyto/6], 'omitnan');
    ER_change = sum([k1_c_ER/6, k2_c_ER/3, k3_c_ER/3, k4_c_ER/6], 'omitnan');
    mito_change = sum([k1_c_mito/6, k2_c_mito/3, k3_c_mito/3, k4_c_mito/6], 'omitnan');
    c_cyto_new = max(0.05, c_cyto + cyto_change);
    c_ER_new = max(0.1, c_ER + ER_change);
    c_mito_new = max(0.05, c_mito + mito_change);

    if ~isreal(c_cyto_new) || ~isreal(c_ER_new) || ~isreal(c_mito_new)
        disp('Complex c_cyto_new')
    end
    %% update E_Ca and h_L
    E_Ca_new = E_Ca_func(c_cyto_new);
    h_L_new = h_L_func(c_cyto_new);
    I_T_type_new = I_T_type_func(V_new, E_Ca_new, m_T_new, h_Tf_new, h_Ts_new);
    I_L_type_new = I_L_type_func(V_new, E_Ca_new, m_L_new, h_L_new);
    I_VGCC_new = I_VGCC_func(I_T_type_new, I_L_type_new);

    %% update all, fill in the record
    IP3 = IP3_new;
    R_particle = R_new;
    PLC_particle = PLC_new;
    DAG = DAG_new;
    act_cPKC_particle = act_cPKC_new;
    AMPA_particle = AMPA_new;
    ATP = ATP_new;
    c_cyto = c_cyto_new;
    c_ER = c_ER_new;
    c_mito = c_mito_new;
    ROS_mito = ROS_mito_new;
    ROS_cyto = ROS_cyto_new;
    glu_neuron = glu_neuron_new;
    glu_cleft = glu_cleft_new;
    glu_ast = glu_ast_new;
    glutamine_ast = glutamine_ast_new;
    
    GsGPCR_particle = GsGPCR_new;
    AC_particle = AC_new;
    cAMP = cAMP_new;
    PKA_particle = PKA_new;
    PDE_particle = PDE_new;

    V = V_new;
    m_T = m_T_new;
    h_Tf = h_Tf_new;
    h_Ts = h_Ts_new;
    m_L = m_L_new;
    E_Ca = E_Ca_new;
    h_L = h_L_new;
    I_T_type = I_T_type_new;
    I_L_type = I_L_type_new;
    I_VGCC = I_VGCC_new;

    IP3_record(ii) = IP3;
    R_record(ii) = mean(R_particle(:,2));
    PLC_record(ii) = mean(PLC_particle(:,2));
    ATP_record(ii) = ATP;
    c_cyto_record(ii) = c_cyto;
    c_ER_record(ii) = c_ER;
    c_mito_record(ii) = c_mito;
    ROS_mito_record(ii) = ROS_mito;
    ROS_cyto_record(ii) = ROS_cyto;
    DAG_record(ii) = DAG;
    act_cPKC_record(ii) = mean(act_cPKC_particle(:,2));
    AMPA_record(ii) = mean(AMPA_particle(:,3));
    glu_neuron_record(ii) = glu_neuron;
    glu_cleft_record(ii) = glu_cleft;
    glu_ast_record(ii) = glu_ast;
    glutamine_ast_record(ii) = glutamine_ast;

    GsGPCR_record(ii) = mean(GsGPCR_new(:,2));
    AC_record(ii) = mean(AC_new(:,2));
    cAMP_record(ii) = cAMP_new;
    PKA_record(ii) = mean(PKA_new(:,2));
    PDE_record(ii) = mean(PDE_new(:,2));

    V_record(ii) = V;
    I_VGCC_record(ii) = I_VGCC;

    FER_record(ii) = b_ER * FER;
    Fserca_record(ii) = b_serca * Fserca;
    F_leakage_record(ii) = b_leakage * F_leakage;
    Fmito_record(ii) = b_Fmito * Fmito;
    Smito_record(ii) = b_Smito * Smito;
    F_in_record(ii) = b_in * F_in;
    F_out_record(ii) = b_out * F_out;
    c_cyto_eff_record(ii) = c_cyto_effect;
    ATP_eff_record(ii) = ATP_effect;
    IP3_eff_record(ii) = IP3_effect(IP3, ATP);
end


chemical_table = table(run_seq', NE_seq', DA_seq', IP3_record', ATP_record', c_cyto_record', c_ER_record', c_mito_record', ROS_mito_record', ROS_cyto_record', R_record', PLC_record', ...
    act_cPKC_record', DAG_record', glu_neuron_record', glu_cleft_record', glu_ast_record', glutamine_ast_record', ...
    GsGPCR_record', AC_record', cAMP_record', PKA_record', PDE_record', AMPA_record', ...
    'VariableNames', {'Running', 'NE', 'DA', 'IP3', 'ATP', 'cytosolic_Ca', 'ER_Ca', 'mito_Ca', 'ROS_mito', 'ROS_cyto', 'Active receptors', 'PLC', ...
    'cPKC*', 'DAG', 'Glu_neuron', 'Glu_cleft', 'Glu_ast', 'Glutamine_ast', 'GsGPCR', 'AC', 'cAMP', 'PKA', 'PDE', 'AMPA'});

flux_table = table(FER_record', Fserca_record', F_leakage_record', Fmito_record', Smito_record', c_cyto_eff_record', ATP_eff_record', IP3_eff_record', F_in_record', F_out_record', V_record', I_VGCC_record', ...
    'VariableNames', {'FER', 'Fserca', 'F_leakage', 'Fmito', 'Smito', 'c_cyto_effect', 'ATP_effect', 'IP3_effect', 'F_in', 'F_out', 'Mem_potential', 'I_VGCC'});
end
