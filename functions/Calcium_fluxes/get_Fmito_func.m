%% modelling the calcium flux from mitochondria to cytosol
% assume the opening of mPTP depends on ROS_mito
% assume the calcium efflux from mito depends on mPTP and c_mito & c_cyto
function Fmito_func = get_Fmito_func(K_fmito, b_mPTP, c_mito, c_cyto, ROS_threshold)
if c_mito > c_cyto
    Fmito_func = @(c_cyto, c_mito, ROS_mito) (c_mito - c_cyto) * mPTP_effect(ROS_mito, c_mito, K_fmito, b_mPTP, ROS_threshold);
else
    Fmito_func = @(c_cyto, c_mito, ROS_mito) 0;
end
end

function mPTP = mPTP_effect(ROS_mito, c_mito, K_fmito, b_mPTP, ROS_threshold)
if ROS_mito >= ROS_threshold || c_mito > 0.5
    % mPTP = b_mPTP * Hill_func(ROS_mito, 1, K_fmito);
    mPTP = 1;
else
    mPTP = 0;
end
end