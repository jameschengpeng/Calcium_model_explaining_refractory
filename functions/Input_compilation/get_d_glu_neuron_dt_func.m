%% model the dynamics of neuronal glutamate
% generation: behavior/stimulation, i.e., the glu_stim
% degradation: exocytosis
function d_glu_neuron_dt_func = get_d_glu_neuron_dt_func(theta, glu_stim_func)
d_glu_neuron_dt_func = @(glu_neuron, glu_cleft_record, t, fps_upsampled) get_time_deri_glu_neuron(glu_neuron, glu_cleft_record, t, fps_upsampled, theta, glu_stim_func);
end

function time_deri_glu_neuron = get_time_deri_glu_neuron(glu_neuron, glu_cleft_record, t, fps_upsampled, theta, glu_stim_func)
b_stim2glu_neuron = theta('b_stim2glu_neuron');
b_glu_neuron2cleft = theta('b_glu_neuron2cleft');
b_glu_cleft2neuron = theta('b_glu_cleft2neuron');
tau_cleft2neuron = theta('tau_cleft2neuron');
if t - tau_cleft2neuron > 0
    idx = ceil(fps_upsampled * (t - tau_cleft2neuron));
    time_deri_glu_neuron = b_stim2glu_neuron * glu_stim_func(t) - b_glu_neuron2cleft * glu_neuron + b_glu_cleft2neuron * glu_cleft_record(idx);
else
    idx = 1;
    time_deri_glu_neuron = b_stim2glu_neuron * glu_stim_func(t) - b_glu_neuron2cleft * glu_neuron + b_glu_cleft2neuron * glu_cleft_record(idx);
end
end