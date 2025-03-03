%% compute the time lag between two seqs, i.e., delay time that yield largest cross correlation
function time_lag_in_second = get_time_lag_in_second(seq1, seq2, fps)
seq1 = reshape(seq1, [], 1);
seq2 = reshape(seq2, [], 1);
delay_time = -10:0.1:10;
shifted_indices = round(delay_time .* fps);
cross_correlation = zeros(size(delay_time));
for ii = 1:length(shifted_indices)
    ss = shifted_indices(ii);
    shifted_seq2 = circshift(seq2, ss);
    if ss >= 0
        cross_correlation(ii) = corr(seq1(ss+1:end), shifted_seq2(ss+1:end));
    else
        cross_correlation(ii) = corr(seq1(1:end+ss), shifted_seq2(1:end+ss));
    end
end
[~, max_ind] = max(cross_correlation);
time_lag_in_second = delay_time(max_ind);
end