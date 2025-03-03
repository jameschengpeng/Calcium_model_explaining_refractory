%% get the dynamics of AMPA receptors
% Saftenku, E. E. "Modeling of slow glutamate diffusion and AMPA receptor activation in the cerebellar glomerulus." Journal of theoretical biology 234.3 (2005): 363-382.

function d_AMPA_particle_dt_func = get_d_AMPA_particle_dt_func(theta)
d_AMPA_particle_dt_func = @(AMPA_particle, glu_cleft) AMPA_particle * get_transition_matrix(glu_cleft, theta);
end

function transition_matrix = get_transition_matrix(glu_cleft, theta)
K_AMPA_Glu = theta('K_AMPA_Glu');
f = glu_cleft^2 / (glu_cleft + K_AMPA_Glu)^2;
ao = 25.39;
ad = 5.11;
ac = 4;
ar = 0.065;
% three states: desensitized, close, open
transition_matrix = [-ar ar 0; ad*f (-ad*f-ao*f) ao*f; 0 ac -ac];
end