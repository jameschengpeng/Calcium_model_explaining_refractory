%% modelling the proportion of active PLC-beta
% PLC is an n*2 matrix, representing the probability distribution of each
% PLC molecule to be either inactive or active
% each PLC molecule is modelled by continuous time Markov chain
function dPLC_particle_dt_func = get_dPLC_particle_dt_func(theta)
% PLC is activated by active receptors, and deactivated by cPKC*
dPLC_particle_dt_func = @(act_cPKC_particle, PLC_particle, R_particle) PLC_particle * get_trans_matrix(act_cPKC_particle, R_particle, theta);
end

% compute the transition probability matrix
function trans_matrix = get_trans_matrix(act_cPKC_particle, R_particle, theta)
b_PLC_produce = theta('b_PLC_produce');
b_PLC_degrade = theta('b_PLC_degrade');
act_cPKC = mean(act_cPKC_particle(:,2));
R = mean(R_particle(:,2));
other_degrade = 0.2;
activate = b_PLC_produce * R; % activated by receptors
deactivate = b_PLC_degrade * act_cPKC + other_degrade; % deactivated by cPKC
trans_matrix = [-activate activate; deactivate -deactivate];
end