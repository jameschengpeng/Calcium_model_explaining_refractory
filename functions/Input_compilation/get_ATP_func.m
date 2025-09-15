%% modelling of ATP synthesis
% return a function 
function ATP_func = get_ATP_func(NE_func, glu_ast, theta)
ATP_func = @(t) theta('b_ATP_NE') * NE_func(t) + theta('b_ATP_glu') * glu_ast;
end