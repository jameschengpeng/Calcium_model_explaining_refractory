%% modelling the production of ROS inside mitochondria
% in: assume the ROS production only depends on glutamate, but there is a delay
% term ROS_delay
% out: assume the ROS is exhausted only via mPTP
function dROS_mito_dt_func = get_dROS_mito_dt_func(theta, fps_upsampled)
ROS_delay = theta('ROS_delay');
b_ROS2glu = theta('b_ROS2glu');
ROS_threshold = theta('ROS_threshold'); % maximum ROS inside mitochondria
b_ROS_flow = theta('b_ROS_flow');

dROS_mito_dt_func = @(ROS_mito, ROS_cyto, c_mito, glu_ast_record, t) b_ROS2glu * ROS_glu(glu_ast_record, ROS_delay, t, fps_upsampled) - b_ROS_flow * ROS_flow(ROS_mito, ROS_cyto, c_mito, ROS_threshold);
end

function dROSdt = ROS_glu(glu_ast_record, ROS_delay, t, fps_upsampled)
if t <= ROS_delay
    dROSdt = 0;
else
    idx = ceil((t-ROS_delay)*fps_upsampled);
    dROSdt = glu_ast_record(idx) - 0.01;
end
end


function flow = ROS_flow(ROS_mito, ROS_cyto, c_mito, ROS_threshold)
if ROS_mito >= ROS_threshold || c_mito > 0.8
    flow = ROS_mito - ROS_cyto;
else
    flow = 0;
end
end