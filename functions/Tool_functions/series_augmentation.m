%% augmenting time series. Augment seq to same length of ref_seq
function aug_seq = series_augmentation(seq, ref_seq)
old_time_course = 0:(length(seq)-1);
new_time_course = 0:(length(ref_seq)-1);
new_time_course = new_time_course.*(length(seq)/length(ref_seq));
new_time_course = min(new_time_course, max(old_time_course));
seq_func = @(t) interp1(old_time_course, seq, t, 'linear');
aug_seq = arrayfun(seq_func, new_time_course);
end