%% compute the dynamics of Gi-GPCR, which is mainly used for capturing gating effect of NE on dopamine
% in mouse astrocytes, GiGPCR contains alpha-2 adrenergic receptors (Adra2a
% gene highly expressed). Though D2-like dopamine receptors are also
% GiGPCR, the DRD2 gene is lowly expressed
function dGiGPCR_particle_dt_func = get_dGiGPCR_particle_dt_func(theta, NE_func, NE_spont, t0, glu_cleft, DA_func)
NE = NE_func(t0);
dGiGPCR_particle_dt_func = @(GiGPCR_particle, PKA_particle, cAMP, t) GiGPCR_particle * get_transition_matrix(theta, NE, NE_spont, NE_func, t, PKA_particle, cAMP, glu_cleft, DA_func);
end

function transition_matrix = get_transition_matrix(theta, NE, NE_spont, NE_func, t, PKA_particle, cAMP, glu_cleft, DA_func)
b_GiGPCR_active = theta('b_GiGPCR_active');
NE_thre = theta('NE_thre');
b_GiGPCR_deactive = theta('b_GiGPCR_deactive');
PKA = mean(PKA_particle(:,2));
if NE > NE_spont + NE_thre
    activate = b_GiGPCR_active * (NE_func(t)-NE_thre + max(0, glu_cleft-0.1) + DA_func(t));
else
    activate = 0;
end
deactivate = b_GiGPCR_deactive * (PKA + cAMP);
transition_matrix = [-activate activate; deactivate -deactivate];
end