% temp_resample_rate = fps_org/fps_upsampled
function [A, B, C, D, E] = get_trial_timing(TT, temp_resample_rate, run_seq)
A = []; B = []; C = []; D = []; E = [];
for ii = 1:size(TT,1)
    trial_type = TT{ii, "trial_type"};
    timing_idx = TT{ii, "end"};
    timing_idx = min(length(run_seq), floor(timing_idx/temp_resample_rate));
    if trial_type == 1
        A = [A timing_idx];
    elseif trial_type == 2
        B = [B timing_idx];
    elseif trial_type == 3
        C = [C timing_idx];
    elseif trial_type == 4
        D = [D timing_idx];
    else
        E = [E timing_idx];
    end
end
end