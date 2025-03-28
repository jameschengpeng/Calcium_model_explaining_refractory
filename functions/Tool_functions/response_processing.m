% given the raw response sequence, i.e., regional average from the image
% data, identify the change of baseline, estimate noise, and compute dF
% resp_fps means the fps for the resp_seq
function [filtered_resp, baseline] = response_processing(resp_seq, resp_fps, cutting_time)
period_length = floor(resp_fps * cutting_time);
total_periods = ceil(length(resp_seq)/period_length);
min_values = zeros(1, total_periods);
for ii = 1:total_periods
    start_idx = (ii-1)*period_length+1;
    end_idx = min(length(resp_seq), ii*period_length);
    subseq = resp_seq(start_idx:end_idx);
    min_values(ii) = get_min_value(subseq);
end
min_values = [min_values(1) min_values];
% linearly interpolate these min values
time_stamps = 0:(length(min_values)-2);
time_stamps = [time_stamps length(resp_seq)/period_length]; % order of each cutting point
time_stamps = time_stamps.*cutting_time; % actual time (s) of each cutting point
min_value_function = @(t) interp1(time_stamps, min_values, t, 'linear');
time_course = 0:(length(resp_seq)-1);
time_course = time_course./resp_fps;
baseline = arrayfun(min_value_function, time_course);
filtered_resp = resp_seq - baseline;
end

function corrected_baseline = get_min_value(seq)
diff_seq = seq(2:end) - seq(1:end-1); % elements of diff_seq are iid. N(0,2*sigma^2)
diff_seq_sq = diff_seq.^2;
sigma_sq = 1/2 * mean(diff_seq_sq);
samples = sqrt(sigma_sq) * randn(length(diff_seq), 1);
sample_min = min(samples);
corrected_baseline = min(seq) - sample_min;
end