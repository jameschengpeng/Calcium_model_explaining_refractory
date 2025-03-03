%% functional form of dV_dt, the membrane potential time derivative
function dV_dt_func = get_dV_dt_func(VGCC_const)
Cm = VGCC_const('Cm');
% dV_dt_func = @(I_VGCC) -1/Cm * I_VGCC;
dV_dt_func = @(I_VGCC) 0;
end