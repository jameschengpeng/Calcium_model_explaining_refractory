%% modelling for PKA
% activated by cAMP. constant degradation
function dPKA_particle_dt_func = get_dPKA_particle_dt_func(theta)
dPKA_particle_dt_func = @(PKA_particle, cAMP) PKA_particle * get_trans_matrix(cAMP, theta);
end

function trans_matrix = get_trans_matrix(cAMP, theta)
b_PKA_produce = theta('b_PKA_produce');
b_PKA_degrade = theta('b_PKA_degrade');
activate = b_PKA_produce * cAMP;
deactivate = b_PKA_degrade;
trans_matrix = [-activate activate; deactivate -deactivate];
end