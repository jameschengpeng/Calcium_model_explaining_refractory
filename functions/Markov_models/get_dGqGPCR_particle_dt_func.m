%% model the dynamics of the percentage of activated Gq-GPCRs
% alpha1-adrenergic receptors, mGluR5, and D1-D2 heteromers are Gq-coupled
function dGqGPCR_particle_dt_func = get_dGqGPCR_particle_dt_func(theta, other_settings, NE_func, DA_func, glu_cleft)
Hill_coeff = 5;
dGqGPCR_particle_dt_func = @(act_cPKC_particle, GqGPCR_particle, GiGPCR_particle, GsGPCR_particle, t) GqGPCR_particle * ...
    get_transition_matrix(act_cPKC_particle, GiGPCR_particle, GsGPCR_particle, ... 
    t, NE_func, DA_func, glu_cleft, Hill_coeff, theta, other_settings);
end

function transition_matrix = get_transition_matrix(act_cPKC_particle, GiGPCR_particle, GsGPCR_particle, t, NE_func, DA_func, ...
    glu_cleft, Hill_coeff, theta, other_settings)
b_GqGPCR_active = theta('b_GqGPCR_active');
b_GqGPCR_deactive = theta('b_GqGPCR_deactive');
NE_thre = theta('NE_thre');
cPKC = mean(act_cPKC_particle(:,2));
GiGPCR = mean(GiGPCR_particle(:,2));
GsGPCR = mean(GsGPCR_particle(:,2));

% Heteromers formed by D1 and D2-like DA receptors are Gq-coupled
% D1-like are Gs-coupled while D2-like are Gi-coupled
DA_coeff = 3 * min(GsGPCR, GiGPCR);

% Gq-GPCR can be activated by NE, DA, or glutamate
activate = b_GqGPCR_active * (max(0, NE_func(t)-NE_thre) + DA_coeff * DA_func(t) + 0.2 * max(0, glu_cleft-0.1));
deactivate = Hill_coeff * cPKC + b_GqGPCR_deactive;

if ~other_settings('Gq_GPCR') % block Gq-GPCR by iBARK
    activate = 0;
end

transition_matrix = [-activate activate; deactivate -deactivate];
end



