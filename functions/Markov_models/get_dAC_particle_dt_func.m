%% modelling the adenylyl cyclase
% model it by particle, and assume that AC is activated by beta-adrenergic
% receptors
% each AC molecule is modelled by a continuous time Markov chain
function dAC_particle_dt_func = get_dAC_particle_dt_func(theta)
dAC_particle_dt_func = @(PKA_particle, AC_particle, GsGPCR_particle) AC_particle * get_trans_matrix(PKA_particle, GsGPCR_particle, theta);
end

% activated by beta-adrenergic receptors (assume proportional to
% alpha-adrenergic recetor
% deactivated by PKA
function trans_matrix = get_trans_matrix(PKA_particle, GsGPCR_particle, theta)
b_AC_produce = theta('b_AC_produce');
b_AC_degrade = theta('b_AC_degrade');
PKA = mean(PKA_particle(:,2));
GsGPCR = mean(GsGPCR_particle(:,2));
activate = b_AC_produce * Hill_func(GsGPCR, 1, 0.2);
if PKA > 0
    deactivate = b_AC_degrade * Hill_func(PKA, 2, 0.3);
else
    deactivate = 0.1;
end
trans_matrix = [-activate activate; deactivate -deactivate];
end