%% functional form of E_Ca, the Nernst potential of Ca
function E_Ca_func = get_E_Ca_func(VGCC_const)
R = VGCC_const('R');
T = VGCC_const('T');
z = VGCC_const('z');
F = VGCC_const('F');
c_ecs = VGCC_const('c_ecs');

% the following line is from the code provided by the original paper
% Eca(n)=59.5/2*log10(Cout(n)/Cin(n)); %%% Eca is Nernst Potential of Ca, 59.5=R*T/F*1000*ln10, unit as mV

E_Ca_func = @(c_cyto)  R*T / (z*F)*1000*log(10) * log10(c_ecs/c_cyto);
end