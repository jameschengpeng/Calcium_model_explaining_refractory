%% model the dynamics of the percentage of activated Gs-GPCRs
% alpha2-adrenergic receptors and D1-like dopamine receptors are Gs-coupled
function dGsGPCR_particle_dt_func = get_dGsGPCR_particle_dt_func(theta, NE_func, NE_spont, DA_func, DA_spont, t0)
NE = NE_func(t0); % t0 is the current time, it should be distinguishable with the general time t
DA = DA_func(t0);
PKA_coeff = 2;
dGsGPCR_particle_dt_func = @(PKA_particle, GsGPCR_particle, t) GsGPCR_particle * get_transition_matrix(PKA_particle, t, NE, NE_func, NE_spont, DA, DA_func, DA_spont, ...
    PKA_coeff, theta);
end

function transition_matrix = get_transition_matrix(PKA_particle, t, NE, NE_func, NE_spont, DA, DA_func, DA_spont, ...
    PKA_coeff, theta)
b_GsGPCR_active = theta('b_GsGPCR_active');
b_GsGPCR_deactive = theta('b_GsGPCR_deactive');
NE_thre = theta('NE_thre');
PKA = mean(PKA_particle(:,2));

% GsGPCR is can be activated by NE or DA, but not by glutamate
if NE > NE_spont + NE_thre
    if DA > DA_spont
        activate = b_GsGPCR_active * (NE_func(t)-NE_thre + 0.6 * DA_func(t));
        deactivate = PKA_coeff * PKA + b_GsGPCR_deactive;
    else
        activate = b_GsGPCR_active * (NE_func(t)-NE_thre);
        deactivate = PKA_coeff * PKA + b_GsGPCR_deactive;
    end
else
    activate = 0;
    deactivate = PKA_coeff * PKA + b_GsGPCR_deactive;
end

transition_matrix = [-activate activate; deactivate -deactivate];
end