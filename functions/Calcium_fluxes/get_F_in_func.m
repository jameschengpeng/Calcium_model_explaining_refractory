%% model the Ca2+ influx via the plasma membrane
% reference the following paper:
% Dupont, Genevieve, and Albert Goldbeter. 
% "One-pool model for Ca2+ oscillations involving Ca2+ and inositol 1, 4, 5-trisphosphate as co-agonists for Ca2+ release." 
% Cell calcium 14.4 (1993): 311-322.

% however, in this article, the authors didn't specify the channel they
% are modelling, so, we turned to VGCC modelling

function F_in_func = get_F_in_func(VGCC_const)
% F_in_func = @(IP3) 0.025 + 0.2 * Hill_func(IP3, 2, 1);
F_in_func = @(I_VGCC, PKA_particle, AMPA_particle) get_J_VGCC(VGCC_const, I_VGCC, PKA_particle) + 0.1 * mean(AMPA_particle(:, 3)) + 0.01; % i.e., returned J_VGCC_func = @(I_VGCC, PKA_particle)
end