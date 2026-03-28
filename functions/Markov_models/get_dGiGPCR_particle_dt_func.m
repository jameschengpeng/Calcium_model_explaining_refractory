%% compute the dynamics of Gi-GPCR, which is mainly used for capturing gating effect of NE on dopamine
% in mouse astrocytes, GiGPCR contains alpha-2 adrenergic receptors (Adra2a
% gene highly expressed). Though D2-like dopamine receptors are also
% GiGPCR, the DRD2 gene is lowly expressed
function dGiGPCR_particle_dt_func = get_dGiGPCR_particle_dt_func(theta, NE_func, NE_spont, DA_func)
dGiGPCR_particle_dt_func = @(GiGPCR_particle, PKA_particle, cAMP, t, glu_cleft) GiGPCR_particle * get_transition_matrix(theta, NE_spont, NE_func, t, PKA_particle, cAMP, glu_cleft, DA_func);
end

function transition_matrix = get_transition_matrix(theta, NE_spont, NE_func, t, PKA_particle, cAMP, glu_cleft, DA_func)
b_GiGPCR_active = theta('b_GiGPCR_active');
NE_thre = theta('NE_thre');
b_GiGPCR_deactive_PKA = theta('b_GiGPCR_deactive_PKA');
b_GiGPCR_deactive_cAMP = theta('b_GiGPCR_deactive_cAMP');
PKA = mean(PKA_particle(:,2));
NE = NE_func(t);
if NE > NE_spont + NE_thre
    activate = b_GiGPCR_active * (theta('NE_Gi') * (NE-NE_thre) + theta('glu_Gi') * max(0, glu_cleft-0.1) + theta('DA_Gi') * DA_func(t));
else
    activate = 0;
end
deactivate = b_GiGPCR_deactive_PKA * PKA + b_GiGPCR_deactive_cAMP * cAMP;
transition_matrix = [-activate activate; deactivate -deactivate];
end