%% functional form of I_T_type
function I_T_type_func = get_I_T_type_func(VGCC_const)
gT_bar = VGCC_const('gT_bar');
I_T_type_func = @(V, E_Ca, m_T, h_Tf, h_Ts) gT_bar * m_T * (h_Tf + 0.04*h_Ts) * (V - E_Ca);
end