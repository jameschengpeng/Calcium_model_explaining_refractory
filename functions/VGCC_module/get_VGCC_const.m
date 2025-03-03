%% function to get the constants in VGCC module
% return a dictionary
function VGCC_const = get_VGCC_const()
VGCC_const = containers.Map;
VGCC_const('R') = 8.31;
VGCC_const('T') = 300;
VGCC_const('z') = 2;
VGCC_const('F') = 96485;
VGCC_const('gT_bar') = 0.06;
VGCC_const('gL_bar') = 3.5;
VGCC_const('V_ast') = 5.233*1e-13; % volume of astrocyte soma, 5 micron radius sphere, in liter
VGCC_const('c_ecs') = 2000; % extracellular calcium concentration in um
VGCC_const('Cm') = 5.3; % membrane capacitance in pF
end