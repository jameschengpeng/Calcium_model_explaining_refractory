%% modelling the ROS in the cytosol
% opening rate of mPTP is proportional to ROS_mito
% in: ROS exits mito via mPTP and enters cytosol
% out: ROS is removed by SOD enzymes in cytosol
function dROS_cyto_dt_func = get_dROS_cyto_dt_func(theta, cyto_mito_ratio, fps_upsampled)
ROS_deletion = theta('ROS_deletion'); % rate of ROS deletion
ROS_threshold = theta('ROS_threshold'); % maximum ROS inside mitochondria
b_ROS_flow = theta('b_ROS_flow');

F_ROS_func = get_F_ROS_func(ROS_threshold, fps_upsampled);
dROS_cyto_dt_func = @(ROS_mito, ROS_cyto, c_mito) b_ROS_flow * F_ROS_func(ROS_mito, ROS_cyto, c_mito)/cyto_mito_ratio - ROS_deletion * ROS_cyto;
end

function F_ROS_func = get_F_ROS_func(ROS_threshold, fps_upsampled)
F_ROS_func = @(ROS_mito, ROS_cyto, c_mito) ROS_flow(ROS_mito, ROS_cyto, ROS_threshold, c_mito, fps_upsampled);
end

function flow = ROS_flow(ROS_mito, ROS_cyto, ROS_threshold, c_mito, fps_upsampled)
sROS = smoothed_step_function(ROS_mito, ROS_threshold, fps_upsampled); % smoothed step function for ROS
sCa  = smoothed_step_function(c_mito, 0.8, fps_upsampled); % smoothed step function for c_mito
s_joint = 1 - (1-sROS)*(1-sCa); % their joint smoothed step function

flow = s_joint * (ROS_mito - ROS_cyto);
end