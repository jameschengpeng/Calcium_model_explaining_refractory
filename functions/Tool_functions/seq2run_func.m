%% given the run_seq and fps information, return the continuous function of run_seq
function run_func = seq2run_func(run_seq, fps_org)
time_course = (0:(length(run_seq)-1))/fps_org;
run_func = @(t) interp1(time_course, run_seq, t, 'linear', 'extrap');
end