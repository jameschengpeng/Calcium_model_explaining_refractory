%% a step function to prevent oscillations near the threshold
function S = smoothed_step_function(chemical, threshold, fps)
dt = 1/fps;
S = 0.5 * (1 + tanh((chemical-threshold)/(dt+eps)));
end