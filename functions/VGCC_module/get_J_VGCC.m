%% functional form of J_VGCC
function J_VGCC = get_J_VGCC(VGCC_const, I_VGCC, PKA_particle)
z = VGCC_const('z');
F = VGCC_const('F');
V_ast = VGCC_const('V_ast');
% this line was from the code provided in supplementary of original paper
% Vin(n)=-in(n)*1e4/(5.23)/(2*F);   %%% Converting the current to the flux with unit as uM/s
J_VGCC = -I_VGCC/(z * F * V_ast * 1e9) * (1 + Hill_func(mean(PKA_particle(:,2)), 1, 0.2));
end