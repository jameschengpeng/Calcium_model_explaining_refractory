%% processing the calcium response sequence
% normalize the sequence to be in range [0,1], where 1 is the maximum value
% of calcium response among ALL trials
function resp_scaled = scale_resp(resp_seq, new_min, new_max)
    if nargin < 3
        new_max = 1;
    end
    if nargin < 2
        new_min = 0.05;
    end

    max_val = max(resp_seq);
    min_val = min(resp_seq);
    resp_scaled = ((resp_seq - min_val) ./ (max_val - min_val)) * (new_max - new_min) + new_min;
end
