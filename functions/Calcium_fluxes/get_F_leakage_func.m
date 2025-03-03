%% modelling the natural leakage of Ca2+ from ER to cytosol
function F_leakage_func = get_F_leakage_func()
% F_leakage_func = @(c_cyto, c_ER) Hill_func(max(0, c_ER-c_cyto), 2, 10);
F_leakage_func = @(c_cyto, c_ER) c_ER - c_cyto;
end