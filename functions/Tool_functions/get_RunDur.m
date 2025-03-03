%% given the velocity record (n_trials * length of running), find the largest running duration
function curr_RunDur = get_RunDur(vv, run_threshold, fps)
[n, ~] = size(vv);
curr_RunDur = zeros(n,1);
for ii = 1:n
    run_seq = vv(ii,:);
    run_seq_filt = imgaussfilt(run_seq, 10); % apply gaussian filtering to obtain lower freq components
    running_periods = run_seq_filt >= run_threshold;
    start_indices = find(diff([0, running_periods]) == 1); % start points of running period
    end_indices = find(diff([running_periods, 0]) == -1); % end points of running period
    run_dur = end_indices - start_indices;
    curr_RunDur(ii) = sum(run_dur)/fps;
end
end