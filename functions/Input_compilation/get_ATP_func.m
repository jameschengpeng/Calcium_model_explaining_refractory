%% modelling of ATP synthesis
% return a function 
function ATP_func = get_ATP_func(NE_func, glu_ast)
ATP_func = @(t) (0.7 * NE_func(t) + 0.3 * glu_ast) * 7;
end