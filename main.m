%% verification on entire recording session
clear;
clc;

% Get the operating system type
arch_type = computer('arch');

% Check operating system and set the base path accordingly
if strcmp(arch_type, 'win64')
    % Windows path format
    base_path = 'C:\Users\james\CBIL\Astrocyte\scalable_calcium_model_prev';
    data_path = 'C:\Users\james\CBIL\Astrocyte\AQUA_processed_data';
elseif strcmp(arch_type, 'glnxa64') 
    % Linux path format
    base_path = '/work/Jamespeng/Astrocyte/scalable_calcium_model_prev';
    data_path = '/work/Jamespeng/Astrocyte/AQUA_processed_data'; 
end

% Add all subdirectories to the MATLAB path
addpath(genpath(base_path));

%%
num = '008';
location = 'northeast';
mouse_name = 'Joey';

loc = 'RecLoc2';
filename = strcat('TAR11R', num, '_', location, '.mat');
aqua_data_path = fullfile(data_path, mouse_name, loc, filename);
fps_upsampled = 50;
fps_org = 30;
[resp_seq, run_seq, stim_onset, reward_onset, TT] = behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled);

%%
init_c_cyto = 0.05; % 50 nM
init_IP3 = 0.25;
init_c_ER = 5; % microM
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
dt = 1/fps_org;
glu_spont = 0.01;
n = 1;
smoothing_factor = 5; % smoothing factor for the run_seq

% three configurations
theta = theta_configuration();
other_settings = other_settings_configuration();
VGCC_const = get_VGCC_const();


[chemical_table, flux_table] = run(theta, other_settings, abs(run_seq), NE_spont, running_threshold, fps_org, fps_upsampled, ...
    stim_onset, cyto_ER_ratio, cyto_mito_ratio, init_IP3, DA_spont, reward_onset, ...
    init_c_cyto, init_c_ER, init_c_mito, init_ROS_mito, init_ROS_cyto, IP3_spont, smoothing_factor, TT);

%%
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
R = chemical_table{:, "Active receptors"};
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

%% simplified plot
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

sgtitle(strcat(num, ' ', location, ' region; ', 'Correlation Coef: ', num2str(corr_coef(1,2)), '; Cosine similarity: ', num2str(cosine_similarity)));




%% simplified plot
figure;
x = (1:length(run_seq)) * 1/fps_upsampled;
ax1 = subplot(3,1,1);
[A, B, C, D, E] = get_trial_timing(TT, fps_org/fps_upsampled, run_seq); % original fps_upsampled = fps_org, by RK4, it is 
p1 = plot(x, filtered_resp, 'b'); % 
hold on; % Hold on to plot more on the same figure
% Plot vertical lines for each list and get legend handles
hA = plotVerticalLines(x, A, 'r');
hB = plotVerticalLines(x, B, 'b');
hC = plotVerticalLines(x, C, '[0, 0.5, 0]');
hD = plotVerticalLines(x, D, 'c');
hE = plotVerticalLines(x, E, '[0.6, 0.4, 0.2]');

% Add legend
% Prepare handles and labels for the legend
handles = [hA, hB, hC, hD, hE];
labels = {'Hit (rewarded!)', 'Miss', 'CR (rewarded!)', 'FA', 'Spont'};
validity = [~isempty(hA) ~isempty(hB) ~isempty(hC) ~isempty(hD) ~isempty(hE)];

xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax2 = gca;
ax2.YAxis(1).Color = 'b';
% Add legend using only valid handles and corresponding labels
legend(handles, labels(validity==1), 'Location', 'westoutside');

% Plot Sequence 1 on the left y-axis
ax2 = subplot(3,1,2);
seq1 = run_seq_filt;
seq2 = filtered_resp;
yyaxis left;
h0 = plot(ax2, x, seq1, 'DisplayName', 'Velocity', 'Color', 'r');
% add vertical line to indicate stimulation and reward if they exist
SA_done = false; 
Thre_done = false;
h11 = [];
h12 = [];
for ii = 1:length(stim_onset)
    if ~ SA_done && stim_onset(ii) > 0
        h11 = xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5, 'DisplayName', 'Salient Visual Stimulation');
        SA_done = true;
    elseif ~ Thre_done && stim_onset(ii) < 0
        h12 = xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5, 'DisplayName', 'Threshold Visual Stimulation');
        Thre_done = true;
    elseif stim_onset(ii) > 0
        xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5);
    else
        xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5);
    end
end

for ii = 1:length(reward_onset)
    if ii == 1
        h2 = xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5, 'DisplayName', 'Reward');
    else
        xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5);
    end
end

xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Velocity', 'FontWeight', 'bold', 'Color', 'r');
% Set the color of the left y-axis to match the curve
ax2 = gca;
ax2.YAxis(1).Color = 'r';
yyaxis right;
h3 = plot(ax2, x, seq2, 'DisplayName', 'Calcium response', 'Color', 'b');
ylabel(ax2, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax2.YAxis(2).Color = 'b';
handles = [h0, h11, h12, h2, h3];
labels = {'Velocity', 'Salient Visual Stimulation', 'Threshold Visual Stimulation', 'Reward', 'Calcium response'};
validity = [~isempty(h0) ~isempty(h11) ~isempty(h12) ~isempty(h2) ~isempty(h3)];
ax2 = gca;
ax2.YAxis(1).Color = 'b';
legend(handles, labels(validity==1), 'Location', 'westoutside');

% Plot Sequence 2,3 on the middle y-axis
ax3 = subplot(3,1,3);
seq3 = filtered_resp;
seq4 = c_cyto;
yyaxis left;
plot(ax3, x, seq3, 'DisplayName', 'Calcium response', 'Color', 'blue');
xlabel(ax3, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax3, 'Calcium response', 'FontWeight', 'bold', 'Color', 'blue');
ax3 = gca;
ax3.YAxis(1).Color = 'blue';
yyaxis right;
plot(ax3, x, seq4, 'DisplayName', 'Predicted response', 'Color', 'black');
ylabel(ax3, 'Predicted response', 'FontWeight', 'bold', 'Color', 'black');
ax3.YAxis(2).Color = 'black';
legend(ax3, 'Location', 'westoutside');

% cosine simularity is invariant on scale of two sequences
cosine_similarity = dot(filtered_resp, c_cyto) / (norm(filtered_resp) * norm(c_cyto));
corr_coef = corrcoef(filtered_resp, c_cyto);

sgtitle(strcat(num, ' ', location, ' region; ', 'Correlation Coef: ', num2str(corr_coef(1,2)), '; Cosine similarity: ', num2str(cosine_similarity)));


%% more contents
figure;
x = (1:length(run_seq)) * 1/fps_upsampled;

ax1 = subplot(7,1,1);
[A, B, C, D, E] = get_trial_timing(TT, fps_org/fps_upsampled, run_seq); % original fps_upsampled = fps_org, by RK4, it is 
p1 = plot(x, filtered_resp, 'b'); % 
hold on; % Hold on to plot more on the same figure
% Plot vertical lines for each list and get legend handles
hA = plotVerticalLines(x, A, 'r');
hB = plotVerticalLines(x, B, 'b');
hC = plotVerticalLines(x, C, '[0, 0.5, 0]');
hD = plotVerticalLines(x, D, 'c');
hE = plotVerticalLines(x, E, '[0.6, 0.4, 0.2]');

% Add legend
% Prepare handles and labels for the legend
handles = [hA, hB, hC, hD, hE];
labels = {'Hit (rewarded!)', 'Miss', 'CR (rewarded!)', 'FA', 'Spont'};
validity = [~isempty(hA) ~isempty(hB) ~isempty(hC) ~isempty(hD) ~isempty(hE)];

xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax2 = gca;
ax2.YAxis(1).Color = 'b';
% Add legend using only valid handles and corresponding labels
legend(handles, labels(validity==1), 'Location', 'westoutside');

% Plot Sequence 1 on the left y-axis
ax2 = subplot(7,1,2);
seq1 = run_seq_filt;
seq2 = filtered_resp;
yyaxis left;
h0 = plot(ax2, x, seq1, 'DisplayName', 'Velocity', 'Color', 'r');
% add vertical line to indicate stimulation and reward if they exist
SA_done = false; 
Thre_done = false;
for ii = 1:length(stim_onset)
    if ~ SA_done && stim_onset(ii) > 0
        h11 = xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5, 'DisplayName', 'Salient Visual Stimulation');
        SA_done = true;
    elseif ~ Thre_done && stim_onset(ii) < 0
        h12 = xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5, 'DisplayName', 'Threshold Visual Stimulation');
        Thre_done = true;
    elseif stim_onset(ii) > 0
        xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5);
    else
        xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5);
    end
end

for ii = 1:length(reward_onset)
    if ii == 1
        h2 = xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5, 'DisplayName', 'Reward');
    else
        xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5);
    end
end

xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Velocity', 'FontWeight', 'bold', 'Color', 'r');
% Set the color of the left y-axis to match the curve
ax2 = gca;
ax2.YAxis(1).Color = 'r';
yyaxis right;
h3 = plot(ax2, x, seq2, 'DisplayName', 'Calcium response', 'Color', 'b');
ylabel(ax2, 'Calcium response', 'FontWeight', 'bold', 'Color', 'b');
ax2.YAxis(2).Color = 'b';
legend(ax2, [h0, h11, h12, h2, h3], 'Location', 'westoutside');

% Plot Sequence 2,3 on the middle y-axis
ax3 = subplot(7,1,3);
seq3 = filtered_resp;
seq4 = c_cyto;
yyaxis left;
plot(ax3, x, seq3, 'DisplayName', 'Calcium response', 'Color', 'blue');
xlabel(ax3, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax3, 'Calcium response', 'FontWeight', 'bold', 'Color', 'blue');
ax3 = gca;
ax3.YAxis(1).Color = 'blue';
yyaxis right;
plot(ax3, x, seq4, 'DisplayName', 'Predicted response', 'Color', 'black');
ylabel(ax3, 'Predicted response', 'FontWeight', 'bold', 'Color', 'black');
ax3.YAxis(2).Color = 'black';
legend(ax3, 'Location', 'westoutside');

% Plot Sequence 1,3 on the last y-axis
ax4 = subplot(7,1,4);
seq5 = cPKC;
seq6 = R;
yyaxis left;
plot(ax4, x, seq5, 'DisplayName', 'cPKC', 'Color', "#D95319");
xlabel(ax4, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax4, 'cPKC', 'FontWeight', 'bold', 'Color', "#D95319");
ax4 = gca;
ax4.YAxis(1).Color = "#D95319";
yyaxis right;
plot(ax4, x, seq6, 'DisplayName', 'Active receptors', 'Color', "#7E2F8E");
ylabel(ax4, 'Active receptors', 'FontWeight', 'bold', 'Color', "#7E2F8E");
ax4.YAxis(2).Color = "#7E2F8E";
legend(ax4, 'Location', 'westoutside');
% Add labels and title

ax5 = subplot(7,1,5);
seq7 = NE;
seq8 = PLC;
yyaxis left;
plot(ax5, x, seq7, 'DisplayName', 'NE', 'Color', "#A2142F");
xlabel(ax5, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax5, 'NE', 'FontWeight', 'bold', 'Color', "#A2142F");
ax5 = gca;
ax5.YAxis(1).Color = "#A2142F";
yyaxis right;
plot(ax5, x, seq8, 'DisplayName', 'PLC', 'Color', 'm');
ylabel(ax5, 'PLC', 'FontWeight', 'bold', 'Color', 'm');
ax5.YAxis(2).Color = 'm';
legend(ax5, 'Location', 'westoutside');

ax6 = subplot(7,1,6);
seq9 = ATP;
seq10 = IP3;
yyaxis left;
plot(ax6, x, seq9, 'DisplayName', 'ATP', 'Color', "#EDB120");
xlabel(ax6, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax6, 'ATP', 'FontWeight', 'bold', 'Color', "#EDB120");
ax6 = gca;
ax6.YAxis(1).Color = "#EDB120";
yyaxis right;
plot(ax6, x, seq10, 'DisplayName', 'IP3', 'Color', "#0072BD");
ylabel(ax6, 'IP3', 'FontWeight', 'bold', 'Color', "#0072BD");
ax6.YAxis(2).Color = "#0072BD";
legend(ax6, 'Location', 'westoutside');


ax7 = subplot(7,1,7);
seq11 = glu_ast;
seq12 = DA;
yyaxis left;
plot(ax7, x, seq11, 'DisplayName', 'Glutamate', 'Color', "#7E2F8E");
xlabel(ax7, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax7, 'Glutamate', 'FontWeight', 'bold', 'Color', "#7E2F8E");
ax7 = gca;
ax7.YAxis(1).Color = "#7E2F8E";
yyaxis right;
plot(ax7, x, seq12, 'DisplayName', 'Dopamine', 'Color', "#77AC30");
ylabel(ax7, 'Dopamine', 'FontWeight', 'bold', 'Color', "#77AC30");
ax7.YAxis(2).Color = "#77AC30";
legend(ax7, 'Location', 'westoutside');

% cosine simularity is invariant on scale of two sequences
cosine_similarity = dot(filtered_resp, c_cyto) / (norm(filtered_resp) * norm(c_cyto));
corr_coef = corrcoef(filtered_resp, c_cyto);

sgtitle(strcat(num, ';', location, ' region; ', 'Correlation Coef: ', num2str(corr_coef(1,2)), '; Cosine similarity: ', num2str(cosine_similarity)));

%% plot for calcium, cAMP, and PKA
figure;
ax1 = subplot(2,1,1);
seq1 = c_cyto;
seq2 = cAMP;
yyaxis left;
plot(ax1, x, seq1, 'DisplayName', 'Predicted Calcium', 'Color', 'black');
xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Predicted Calcium', 'FontWeight', 'bold', 'Color', 'black');
ax1 = gca;
ax1.YAxis(1).Color = 'black';
yyaxis right;
plot(ax1, x, seq2, 'DisplayName', 'cAMP', 'Color', 'red');
ylabel(ax1, 'cAMP', 'FontWeight', 'bold', 'Color', 'red');
ax1.YAxis(2).Color = 'red';
legend(ax1, 'Location', 'westoutside');

ax2 = subplot(2,1,2);
seq3 = c_cyto;
seq4 = PKA;
yyaxis left;
plot(ax2, x, seq3, 'DisplayName', 'Predicted Calcium', 'Color', 'black');
xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Predicted Calcium', 'FontWeight', 'bold', 'Color', 'black');
ax2 = gca;
ax2.YAxis(1).Color = 'black';
yyaxis right;
plot(ax2, x, seq4, 'DisplayName', 'PKA', 'Color', 'magenta');
ylabel(ax2, 'PKA', 'FontWeight', 'bold', 'Color', 'magenta');
ax2.YAxis(2).Color = 'magenta';
legend(ax2, 'Location', 'westoutside');


%% plot NTs vs running
figure;
ax1 = subplot(3,1,1);
seq1 = run_seq;
seq2 = NE;
yyaxis left;
h0 = plot(ax1, x, seq1, 'DisplayName', 'Velocity', 'Color', '#77AC30');
% add vertical line to indicate stimulation and reward if they exist
SA_done = false; 
Thre_done = false;
for ii = 1:length(stim_onset)
    if ~ SA_done && stim_onset(ii) > 0
        h11 = xline(ax1, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5, 'DisplayName', 'Salient Visual Stimulation');
        SA_done = true;
    elseif ~ Thre_done && stim_onset(ii) < 0
        h12 = xline(ax1, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5, 'DisplayName', 'Threshold Visual Stimulation');
        Thre_done = true;
    elseif stim_onset(ii) > 0
        xline(ax1, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5);
    else
        xline(ax1, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5);
    end
end

for ii = 1:length(reward_onset)
    if ii == 1
        h2 = xline(ax1, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5, 'DisplayName', 'Reward');
    else
        xline(ax1, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5);
    end
end

xlabel(ax1, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax1, 'Velocity', 'FontWeight', 'bold', 'Color', '#77AC30');
% Set the color of the left y-axis to match the curve
ax1 = gca;
ax1.YAxis(1).Color = '#77AC30';
yyaxis right;
h3 = plot(ax1, x, seq2, 'DisplayName', 'NE', 'Color', 'k');
ylabel(ax1, 'NE', 'FontWeight', 'bold', 'Color', 'k');
ax1.YAxis(2).Color = 'k';
legend(ax1, [h0, h11, h12, h2, h3], 'Location', 'westoutside');

%%%%%%%%%%%%%%%%%%%%%%% Plot Sequence 1 on the left y-axis
ax2 = subplot(3,1,2);
seq3 = run_seq;
seq4 = glu_ast;
yyaxis left;
h0 = plot(ax2, x, seq3, 'DisplayName', 'Velocity', 'Color', '#77AC30');
% add vertical line to indicate stimulation and reward if they exist
SA_done = false; 
Thre_done = false;
h21 = [];
h22 = [];
for ii = 1:length(stim_onset)
    if ~ SA_done && stim_onset(ii) > 0
        h21 = xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5, 'DisplayName', 'Salient Visual Stimulation');
        SA_done = true;
    elseif ~ Thre_done && stim_onset(ii) < 0
        h22 = xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5, 'DisplayName', 'Threshold Visual Stimulation');
        Thre_done = true;
    elseif stim_onset(ii) > 0
        xline(ax2, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5);
    else
        xline(ax2, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5);
    end
end

for ii = 1:length(reward_onset)
    if ii == 1
        h2 = xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5, 'DisplayName', 'Reward');
    else
        xline(ax2, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5);
    end
end

xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Velocity', 'FontWeight', 'bold', 'Color', '#77AC30');
% Set the color of the left y-axis to match the curve
ax2 = gca;
ax2.YAxis(1).Color = '#77AC30';
yyaxis right;
h3 = plot(ax2, x, seq4, 'DisplayName', 'Glutamate', 'Color', '#7E2F8E');
ylabel(ax2, 'Glutamate', 'FontWeight', 'bold', 'Color', '#7E2F8E');
ax2.YAxis(2).Color = '#7E2F8E';
handles = [h0, h21, h22, h2, h3];
labels = {'Velocity', 'Salient Visual Stimulation', 'Threshold Visual Stimulation', 'Reward', 'Calcium response'};
validity = [~isempty(h0) ~isempty(h11) ~isempty(h12) ~isempty(h2) ~isempty(h3)];
ax2 = gca;
ax2.YAxis(1).Color = '#7E2F8E';
legend(handles, labels(validity==1), 'Location', 'westoutside');

% Plot Sequence 2,3 on the middle y-axis
ax3 = subplot(3,1,3);
seq5 = run_seq;
seq6 = DA;
yyaxis left;
h0 = plot(ax3, x, seq5, 'DisplayName', 'Velocity', 'Color', '#77AC30');
% add vertical line to indicate stimulation and reward if they exist
SA_done = false; 
Thre_done = false;
h31 = [];
h32 = [];
for ii = 1:length(stim_onset)
    if ~ SA_done && stim_onset(ii) > 0
        h31 = xline(ax3, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5, 'DisplayName', 'Salient Visual Stimulation');
        SA_done = true;
    elseif ~ Thre_done && stim_onset(ii) < 0
        h32 = xline(ax3, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5, 'DisplayName', 'Threshold Visual Stimulation');
        Thre_done = true;
    elseif stim_onset(ii) > 0
        xline(ax3, stim_onset(ii)/fps_org, 'r', 'LineWidth', 0.5);
    else
        xline(ax3, abs(stim_onset(ii))/fps_org, '--r', 'LineWidth', 0.5);
    end
end

for ii = 1:length(reward_onset)
    if ii == 1
        h2 = xline(ax3, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5, 'DisplayName', 'Reward');
    else
        xline(ax3, reward_onset(ii)/fps_org, 'b', 'LineWidth', 0.5);
    end
end

xlabel(ax3, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax3, 'Velocity', 'FontWeight', 'bold', 'Color', '#77AC30');
% Set the color of the left y-axis to match the curve
ax3 = gca;
ax3.YAxis(1).Color = '#77AC30';
yyaxis right;
h3 = plot(ax3, x, seq6, 'DisplayName', 'DA', 'Color', '#D95319');
ylabel(ax3, 'DA', 'FontWeight', 'bold', 'Color', '#D95319');
ax3.YAxis(2).Color = '#D95319';
handles = [h0, h31, h32, h2, h3];
labels = {'Velocity', 'Salient Visual Stimulation', 'Threshold Visual Stimulation', 'Reward', 'Calcium response'};
validity = [~isempty(h0) ~isempty(h11) ~isempty(h12) ~isempty(h2) ~isempty(h3)];
ax2 = gca;
ax2.YAxis(1).Color = '#D95319';
legend(handles, labels(validity==1), 'Location', 'westoutside');

sgtitle('NE, Glutamate, DA');

%% save result
save_result = true;
if save_result
    savepath = fullfile('C:\Users\james\CBIL\Astrocyte\Temporal_model_eval', mouse_name, loc, strcat('TAR10R', num));
    savefile = fullfile(savepath, strcat('Results_', num2str(corr_coef(1,2)), '.mat'));
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    save(savefile, 'chemical_table', 'flux_table', 'full_vel', 'TT', 'regional_response_temp_down');
end

%%
real_trace = resp_seq';
pred_trace = c_cyto;
sampling_rate = fps_upsampled;
% Compute the FFT of the signals
real_fft = fft(real_trace);
pred_fft = fft(pred_trace);

% Length of the signal
N = length(real_trace);

% Number of unique points
numUniquePts = ceil((N+1)/2);

% Compute the power spectrum
real_power = abs(real_fft(1:numUniquePts)).^2;
pred_power = abs(pred_fft(1:numUniquePts)).^2;

% Normalize the power spectra
real_power = real_power / sum(real_power);
pred_power = pred_power / sum(pred_power);

% Frequency axis
freq = (0:numUniquePts-1) * (sampling_rate / N);

% Take the logarithm of the power values
log_real_power = log10(real_power);
log_pred_power = log10(pred_power);

% Plot the log power spectra
figure;
plot(freq, log_real_power, 'r-', 'DisplayName', 'Real Power Spectrum');
hold on;
plot(freq, log_pred_power, 'b--', 'DisplayName', 'Predicted Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Log Power');
legend;
title('Logarithmic Power Spectrum Comparison');



%%
% Function to plot vertical lines and create a handle for the legend
function h = plotVerticalLines(x, indices, color)
    if isempty(indices)
        h = []; % No handle if the list is empty
    else
        xline(x(indices(1)), 'Color', color, 'LineWidth', 1); % Plot the first line for legend
        h = plot(NaN, NaN, 'Color', color, 'LineWidth', 1); % Invisible line for legend
        for index = indices(2:end) % Start from second element to avoid duplicating first line
            xline(x(index), 'Color', color, 'LineWidth', 1);
        end
    end
end
