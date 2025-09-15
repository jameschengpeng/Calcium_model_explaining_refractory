%% model the dynamics of glutamine
% generation: from glu_ast
% degradation: exocytosis
function d_glutamine_ast_dt_func = get_d_glutamine_ast_dt_func(theta)
b_glu_ast2glutamine = theta('b_glu_ast2glutamine');
b_glu_glutamine2neu = theta('b_glu_glutamine2neu');
d_glutamine_ast_dt_func = @(glu_ast, glutamine_ast) b_glu_ast2glutamine * glu_ast - b_glu_glutamine2neu * glutamine_ast;
end