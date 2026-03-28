%% modelling of ATP synthesis
% return a function 
function ATP_func = get_ATP_func(NE_func, theta)
b_ATP_NE = theta('b_ATP_NE');
b_ATP_glu = theta('b_ATP_glu');
ATP_func = @(t, glu_ast) b_ATP_NE * NE_func(t) + b_ATP_glu * glu_ast;
end