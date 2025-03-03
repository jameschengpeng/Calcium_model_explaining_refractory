%% modelling the formation of IP3, use the time t to do Runge-Kutta approximation
% return a function handle with ONLY IP3 and t as inputs (i.e., in the form
% dy/dx = f(x,y)). In other words, return the partial derivative of IP3
% with respect to time
% note that dIP3_dt depends on IP3 itself
function dIP3_dt_func = get_dIP3_dt_func(theta)
dIP3_dt_func = @(c_cyto, IP3, PLC_particle) get_time_deri_IP3(c_cyto, IP3, PLC_particle, theta);
end

function time_deri_IP3 = get_time_deri_IP3(c_cyto, IP3, PLC_particle, theta)
b_IP3 = theta('b_IP3');
b_IP3_degrade1 = theta('b_IP3_degrade1');
b_IP3_degrade2 = theta('b_IP3_degrade2');
K_IP3_degrade = theta('K_IP3_degrade');
K_IP3_degrade_c_cyto = theta('K_IP3_degrade_c_cyto');
K_IP3_DAG_prod = theta('K_IP3_DAG_prod');
PLC = mean(PLC_particle(:,2));
% the degradation references Eqn 5.29 in the book CGS, note that the
% affinity factor of c_cyto depends on calmodulin, so not necessarily a
% constant
if c_cyto >= 0.05
    % time_deri_IP3 = b_IP3 * R * reverse_Hill_func(PLC, 2, 0.1) - b_IP3_degrade *  Hill_func(c_cyto, 4, K_IP3_degrade_c_cyto/c_cyto) * Hill_func(IP3, 1, K_IP3_degrade);
    time_deri_IP3 = b_IP3 * Hill_func(PLC, 1, K_IP3_DAG_prod) - b_IP3_degrade1 *  Hill_func(c_cyto, 4, K_IP3_degrade_c_cyto/(1+c_cyto)) * Hill_func(IP3, 1, K_IP3_degrade) ...
        - b_IP3_degrade2 * IP3;
else
    time_deri_IP3 = - b_IP3_degrade2 * IP3;
end
end