%% resample seq1 to same length as seq2
function seq1_resampled = resample_seq(seq1, seq2)
n1 = numel(seq1);
n2 = numel(seq2);
xi = linspace(1, n1, n2);
seq1_resampled = interp1(1:n1, seq1, xi, 'cubic');
end