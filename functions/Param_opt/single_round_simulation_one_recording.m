%% perform one round simulation on one recording.
function [chemical_table, flux_table] = single_round_simulation_one_recording(theta, other_settings, run_seq, stim_onset, reward_onset, TT, ...
                precomputed_NE_seq, precomputed_NE_func, precomputed_glu_stim_seq, precomputed_glu_stim_func, precomputed_DA_seq, precomputed_DA_func)
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

[chemical_table, flux_table] = run(theta, other_settings, abs(run_seq), NE_spont, running_threshold, fps_org, fps_upsampled, ...
    stim_onset, cyto_ER_ratio, cyto_mito_ratio, init_IP3, DA_spont, reward_onset, ...
    init_c_cyto, init_c_ER, init_c_mito, init_ROS_mito, init_ROS_cyto, IP3_spont, smoothing_factor, TT, ...
    use_precomputed_input, precomputed_NE_seq, precomputed_NE_func, precomputed_glu_stim_seq, precomputed_glu_stim_func, precomputed_DA_seq, precomputed_DA_func, ...
    use_neuron_derived_glu, include_ATP);
end