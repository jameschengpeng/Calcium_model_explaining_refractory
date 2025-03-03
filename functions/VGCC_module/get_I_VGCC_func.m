%% function to compute I_VGCC, summation of T-type and L-type VGCC
function I_VGCC_func = get_I_VGCC_func()
I_VGCC_func = @(I_T_type, I_L_type) 2 * I_T_type + I_L_type;
end