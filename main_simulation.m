%% add path. NO NEED TO CHANGE ANYTHING HERE.
clear;
clc;

base_path = fileparts(which('main_simulation.m'));
data_path = fullfile(base_path, 'Data');

% Add all subdirectories to the MATLAB path
addpath(genpath(base_path));

%% Specify the mouse name, recording location, and recording number. FEEL FREE TO CHANGE HERE BASED ON THE DATA FOLDER
mouse_name = 'Joey'; % This can be changed. See Data/
loc = 'RecLoc2'; % This can be changed. See Data/
num = '010'; % This can be changed. See Data/
prefix = strcat(mouse_prefix(mouse_name), num);
files = dir(fullfile(data_path, mouse_name, loc, [prefix '*.mat']));
aqua_data_path = fullfile(data_path, mouse_name, loc, files(1).name);
fps_upsampled = 50;
fps_org = 30;
[resp_seq, run_seq, stim_onset, reward_onset, TT] = behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled);

%% Doing the simulation. These initial values are same as what is presented in the paper.
init_c_cyto = 0.05; % 50 nM
init_IP3 = 0.25;
init_c_ER = 50; % microM
init_c_mito = init_c_cyto; % initially, the mitochondrial calcium is same as cytosolic calcium
init_ROS_mito = 0.1;
init_ROS_cyto = 0;

% constants
cyto_ER_ratio = 5;
cyto_mito_ratio = 20;
NE_spont = 0.01;
DA_spont = 0.01;
IP3_spont = init_IP3;
running_threshold = 0.5;
fps_upsampled = 50;
fps_org = 30;
smoothing_factor = 5; % smoothing factor for the run_seq

use_precomputed_input = true; % 
use_neuron_derived_glu = false; 
include_ATP = true;

theta = theta_configuration();
other_settings = other_settings_configuration();

[pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun] = ...
    precompute_inputs(theta, other_settings, aqua_data_path, fps_org, fps_upsampled);


[chemical_table, flux_table] = run(theta, other_settings, abs(run_seq), NE_spont, running_threshold, fps_org, fps_upsampled, ...
    stim_onset, cyto_ER_ratio, cyto_mito_ratio, init_IP3, DA_spont, reward_onset, ...
    init_c_cyto, init_c_ER, init_c_mito, init_ROS_mito, init_ROS_cyto, IP3_spont, smoothing_factor, TT, ...
    use_precomputed_input, pre_NE_seq, pre_NE_fun, pre_glu_seq, pre_glu_fun, pre_DA_seq, pre_DA_fun, ...
    use_neuron_derived_glu, include_ATP);

%% Extract model predictions
FER = flux_table{:,"FER"};
Fserca = flux_table{:,"Fserca"};
F_leakage = flux_table{:,"F_leakage"};
Fmito = flux_table{:,"Fmito"};
Smito = flux_table{:,"Smito"};
F_in = flux_table{:, "F_in"};
F_out = flux_table{:, "F_out"};
mem_potential = flux_table{:, "Mem_potential"};
I_VGCC = flux_table{:, "I_VGCC"};
c_cyto_effect = flux_table{:, "c_cyto_effect"};
ATP_effect = flux_table{:, "ATP_effect"};
IP3_effect = flux_table{:, "IP3_effect"};
NE = chemical_table{:,"NE"};
DA = chemical_table{:,"DA"};
IP3 = chemical_table{:,"IP3"};
ATP = chemical_table{:,"ATP"};
c_cyto = chemical_table{:,"cytosolic_Ca"};
c_ER = chemical_table{:,"ER_Ca"};
c_mito = chemical_table{:, "mito_Ca"};
ros_mito = chemical_table{:, "ROS_mito"};
ros_cyto = chemical_table{:, "ROS_cyto"};
GqGPCR = chemical_table{:, "GqGPCR"};
PLC = chemical_table{:, "PLC"};
cPKC = chemical_table{:, "cPKC*"};
DAG = chemical_table{:, "DAG"};
glu_neuron = chemical_table{:, "Glu_neuron"};
glu_cleft = chemical_table{:, "Glu_cleft"};
glu_ast = chemical_table{:, "Glu_ast"};
glutamine_ast = chemical_table{:, "Glutamine_ast"};
AMPA = chemical_table{:, "AMPA"};


GsGPCR = chemical_table{:, "GsGPCR"};
AC = chemical_table{:, "AC"};
cAMP = chemical_table{:, "cAMP"};
PKA = chemical_table{:, "PKA"};
PDE = chemical_table{:, "PDE"};

%
run_seq = chemical_table{:,'Running'};
run_seq_filt = imgaussfilt(run_seq, 1);


%% Plot the prediction
figure;
x = (1:length(run_seq)) * 1/fps_upsampled;
filtered_resp = imgaussfilt(resp_seq, fps_upsampled);
filtered_resp = scale_resp(filtered_resp);
ax1 = subplot(3,1,1);
seq1 = run_seq_filt;
seq2 = filtered_resp;
yyaxis left;
plot(ax1, x, seq1, 'DisplayName', 'Velocity', 'Color', 'r');
xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Velocity', 'FontWeight', 'bold', 'Color', 'r');
% Set the color of the left y-axis to match the curve
ax1 = gca;
ax1.YAxis(1).Color = 'r';
yyaxis right
plot(x, seq2, 'b'); % 
xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax1 = gca;
ax1.YAxis(2).Color = 'b';


% Plot Sequence 1 on the left y-axis
ax2 = subplot(3,1,2);
seq3 = c_cyto;
seq4 = filtered_resp;
yyaxis left;
plot(ax2, x, seq3, 'DisplayName', 'Predicted response', 'Color', 'k');
xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Predicted response', 'FontWeight', 'bold', 'Color', 'k');
% Set the color of the left y-axis to match the curve
ax2 = gca;
ax2.YAxis(1).Color = 'k';
yyaxis right
plot(x, seq4, 'b'); % 
xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax2 = gca;
ax2.YAxis(2).Color = 'b';

% Plot Sequence 2,3 on the middle y-axis
ax3 = subplot(3,1,3);
seq5 = c_cyto;
seq6 = cPKC;
yyaxis left;
plot(ax3, x, seq5, 'DisplayName', 'Calcium response', 'Color', 'k');
xlabel(ax3, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax3, 'Predicted response', 'FontWeight', 'bold', 'Color', 'k');
ax3 = gca;
ax3.YAxis(1).Color = 'k';
yyaxis right;
plot(ax3, x, seq6, 'DisplayName', 'cPKC', 'Color', '#D95319');
ylabel(ax3, 'cPKC', 'FontWeight', 'bold', 'Color', '#D95319');
ax3.YAxis(2).Color = '#D95319';


% cosine simularity is invariant on scale of two sequences
cosine_similarity = dot(filtered_resp, c_cyto) / (norm(filtered_resp) * norm(c_cyto));
corr_coef = corrcoef(filtered_resp, c_cyto);

sgtitle(strcat('Correlation Coef: ', num2str(corr_coef(1,2)), '; Cosine similarity: ', num2str(cosine_similarity)));

