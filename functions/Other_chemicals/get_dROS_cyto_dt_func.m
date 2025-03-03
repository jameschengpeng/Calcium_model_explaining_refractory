%% modelling the ROS in the cytosol
% opening rate of mPTP is proportional to ROS_mito
% in: ROS exits mito via mPTP and enters cytosol
% out: ROS is removed by SOD enzymes in cytosol
function dROS_cyto_dt_func = get_dROS_cyto_dt_func(theta, cyto_mito_ratio)
ROS_deletion = theta('ROS_deletion'); % rate of ROS deletion
ROS_threshold = theta('ROS_threshold'); % maximum ROS inside mitochondria
b_ROS_flow = theta('b_ROS_flow');

F_ROS_func = get_F_ROS_func(ROS_threshold);
dROS_cyto_dt_func = @(ROS_mito, ROS_cyto) b_ROS_flow * F_ROS_func(ROS_mito, ROS_cyto)/cyto_mito_ratio - ROS_deletion * ROS_cyto;
end

function F_ROS_func = get_F_ROS_func(ROS_threshold)
F_ROS_func = @(ROS_mito, ROS_cyto) ROS_flow(ROS_mito, ROS_cyto, ROS_threshold);
end

function flow = ROS_flow(ROS_mito, ROS_cyto, ROS_threshold)
if ROS_mito >= ROS_threshold
    flow = ROS_mito - ROS_cyto;
else
    flow = 0;
end
end