%% modelling for the proportion of cPKC*, within range 0~1
function d_act_cPKC_particle_dt_func = get_d_act_cPKC_particle_dt_func(theta)
% cPKC* is activated by one Ca2+ and one DAG molecule
d_act_cPKC_particle_dt_func = @(c_cyto, c_cyto_record, DAG, DAG_record, act_cPKC_particle, fps_upsampled, ii) act_cPKC_particle * get_transition_matrix(c_cyto, c_cyto_record, DAG, DAG_record, theta, fps_upsampled, ii, act_cPKC_particle);
end

function transition_matrix = get_transition_matrix(c_cyto, c_cyto_record, DAG, DAG_record, theta, fps_upsampled, ii, act_cPKC_particle)
b_cpkc_produce = theta('b_cpkc_produce');
b_cpkc_degrade = theta('b_cpkc_degrade');
tau_cpkc_degrade = theta('tau_cpkc_degrade');
c_cyto_peak = theta('c_cyto_peak');
cPKC_delay = theta('cPKC_delay');
if cPKC_delay ~= 0
    ref_idx = max(1, ii - round(cPKC_delay * fps_upsampled));
    c_cyto_ref = c_cyto_record(ref_idx);
    DAG_ref = DAG_record(ref_idx);
else
    c_cyto_ref = c_cyto;
    DAG_ref = DAG;
end

Kc_cyto = 0.6;
KDAG = 0.2;
local_start = max(1, ii - round(tau_cpkc_degrade * fps_upsampled));
local_c_cyto_seq = c_cyto_record(local_start:ii);
marker_point = min(max(max(local_c_cyto_seq), 1), c_cyto_peak);
local_closest2marker = find(abs(local_c_cyto_seq-marker_point) <= 0.03, 1, 'last' );
c_cyto_thre_idx = local_start-1 + local_closest2marker;
if (~isempty(c_cyto_thre_idx) && (ii - c_cyto_thre_idx)/fps_upsampled <= tau_cpkc_degrade)
    activate = b_cpkc_produce * Hill_func(c_cyto_ref, 1, Kc_cyto) * Hill_func(DAG_ref, 1, KDAG);
    deactivate = 0.01 * b_cpkc_degrade;
    % deactivate = 0;
else
    activate = b_cpkc_produce * Hill_func(c_cyto_ref, 1, Kc_cyto) * Hill_func(DAG_ref, 1, KDAG);
    % deactivate = b_cpkc_degrade * reverse_Hill_func(c_cyto, 2, 0.6);
    deactivate = b_cpkc_degrade;
end

transition_matrix = [-activate activate; deactivate -deactivate];
end