%% modelling the dynamics of phosphodiesterase (PDE)
% PDE is activated by PKA, constantly degrade
function dPDE_particle_dt_func = get_dPDE_particle_dt_func(theta)
dPDE_particle_dt_func = @(PDE_particle, c_cyto) PDE_particle * get_trans_matrix(c_cyto, PDE_particle, theta);
end

function trans_matrix = get_trans_matrix(c_cyto, PDE_particle, theta)
b_PDE_produce = theta('b_PDE_produce');
b_PDE_degrade = theta('b_PDE_degrade');
PDE = mean(PDE_particle(:,2));
activate = b_PDE_produce * Hill_func(c_cyto, 4, 0.5); % activated by CaM Ca4 complex 
degrade = b_PDE_degrade;
trans_matrix = [-activate activate; degrade -degrade];
end