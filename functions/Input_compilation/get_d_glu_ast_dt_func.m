%% model the glutamate in the astrocytes
% generation: from glu_cleft via the transporters
% degradation: to glutamine
function d_glu_ast_dt_func = get_d_glu_ast_dt_func(theta)
b_glu_cleft2ast = theta('b_glu_cleft2ast');
b_glu_ast2glutamine = theta('b_glu_ast2glutamine');
d_glu_ast_dt_func = @(glu_ast, glu_cleft) b_glu_cleft2ast * glu_cleft - b_glu_ast2glutamine * glu_ast;
end