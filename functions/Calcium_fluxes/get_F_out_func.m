%% modelling the calcium leaving cytosol to extracellular space
% reference to the following paper:
% Ullah, Ghanim, Peter Jung, and Ann H. Cornell-Bell. 
% "Anti-phase calcium oscillations in astrocytes via inositol (1, 4, 5)-trisphosphate regeneration." 
% Cell calcium 39.3 (2006): 197-208.

function F_out_func = get_F_out_func()
% F_out_func = @(c_cyto, act_cPKC_particle) 0.1 * c_cyto + 0.1 * mean(act_cPKC_particle(:,2));
F_out_func = @(c_cyto, act_cPKC_particle) c_cyto * (0.7 + 2 * mean(act_cPKC_particle(:,2)));
end