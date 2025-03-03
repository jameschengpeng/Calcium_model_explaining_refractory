%% functional form of I_L_type
function I_L_type_func = get_I_L_type_func(VGCC_const)
gL_bar = VGCC_const('gL_bar');
I_L_type_func = @(V, E_Ca, m_L, h_L) gL_bar * m_L * h_L * (V - E_Ca);
end