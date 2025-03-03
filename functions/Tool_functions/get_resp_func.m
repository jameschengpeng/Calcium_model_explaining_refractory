%% get the function of response sequence
function resp_func = get_resp_func(resp_seq, fps_org)
time_course = (0:(length(resp_seq)-1))/fps_org;
resp_func = @(t) interp1(time_course, resp_seq, t, 'linear');
end