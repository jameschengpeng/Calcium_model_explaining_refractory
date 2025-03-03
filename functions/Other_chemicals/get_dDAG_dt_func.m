%% model the dynamics of DAG
function dDAG_dt_func = get_dDAG_dt_func(theta)
% synthesis is coupled with IP3, degradation is mediated by DAG kinase
dDAG_dt_func = @(c_cyto, PLC_particle, DAG) get_time_deri_DAG(theta, c_cyto, PLC_particle, DAG);
end

function time_deri_DAG = get_time_deri_DAG(theta, c_cyto, PLC_particle, DAG)
b_IP3 = theta('b_IP3'); % production rate of DAG is the same as IP3
b_DAG_degrade = theta('b_DAG_degrade');
K_IP3_DAG_prod = theta('K_IP3_DAG_prod');
PLC = mean(PLC_particle(:,2));
other_degrade = 2;
if c_cyto > 0.05
    time_deri_DAG = b_IP3 * Hill_func(PLC, 1, K_IP3_DAG_prod) - b_DAG_degrade * Hill_func(c_cyto, 2, 1) * Hill_func(DAG, 2, 1) - other_degrade * DAG;
else
    time_deri_DAG = - other_degrade * DAG;
end
end