%% model the glutamate in the synaptic cleft
% generation: from glu_neuron
% degradation: return to glu_neuron, or enter astrocyte by transporters
function d_glu_cleft_dt_func = get_d_glu_cleft_dt_func(theta)
d_glu_cleft_dt_func = @(glu_cleft_record, glu_cleft, glu_neuron, t, fps_upsampled) get_time_deri_glu_cleft(theta, glu_cleft_record, glu_cleft, glu_neuron, t, fps_upsampled);
end

function time_deri_glu_cleft = get_time_deri_glu_cleft(theta, glu_cleft_record, glu_cleft, glu_neuron, t, fps_upsampled)
b_glu_neuron2cleft = theta('b_glu_neuron2cleft');
b_glu_cleft2neuron = theta('b_glu_cleft2neuron');
b_glu_cleft2ast = theta('b_glu_cleft2ast');
tau_cleft2neuron = 1;
if t - tau_cleft2neuron > 0
    idx = ceil(fps_upsampled * (t - tau_cleft2neuron));
    % exocytosis from neuron; going back to neuron; entering astrocytes
    time_deri_glu_cleft = b_glu_neuron2cleft * glu_neuron - b_glu_cleft2neuron * glu_cleft_record(idx) - b_glu_cleft2ast * glu_cleft;
else
    idx = 1;
    time_deri_glu_cleft = b_glu_neuron2cleft * glu_neuron - b_glu_cleft2neuron * glu_cleft_record(idx) - b_glu_cleft2ast * glu_cleft;
end
end