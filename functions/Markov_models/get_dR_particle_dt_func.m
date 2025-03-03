%% model the dynamics of the percentage of activated Gq-GPCRs
% use R for activated receptor, r for unactivated
function dR_particle_dt_func = get_dR_particle_dt_func(theta, other_settings, NE_func, ~, DA_func, ~, glu_cleft, t0)
NE = NE_func(t0); % t0 is the current time, it should be distinguishable with the general time t
DA = DA_func(t0);
Hill_coeff = 5;
other_degrade = 0.8;
dR_particle_dt_func = @(act_cPKC_particle, R_particle, t) R_particle * get_transition_matrix(act_cPKC_particle, t, NE_func, DA_func, ...
    glu_cleft, Hill_coeff, other_degrade, theta, other_settings);
end

function transition_matrix = get_transition_matrix(act_cPKC_particle, t, NE_func, DA_func, ...
    glu_cleft, Hill_coeff, other_degrade, theta, other_settings)
b_recep_active = theta('b_recep_active');
NE_thre = theta('NE_thre');
cPKC = mean(act_cPKC_particle(:,2));

% Gq-GPCR can be activated by NE, DA, or glutamate
activate = b_recep_active * (max(0, NE_func(t)-NE_thre) + max(0, 0.3 * DA_func(t)) + 0.2 * max(0, glu_cleft-0.1));
deactivate = Hill_coeff * cPKC + other_degrade;

if ~other_settings('Gq_GPCR') % block Gq-GPCR by iBARK
    activate = 0;
end

transition_matrix = [-activate activate; deactivate -deactivate];
end



