%% evaluate real and predicted signals (concatenated) over multiple metrics
function metric_vector = multi_metric_eval(simulation_result, real_result, save_path)
fields = fieldnames(simulation_result);
corr = zeros(length(fields), 1);
cos = zeros(length(fields), 1);
rmse = zeros(length(fields), 1);
rmse_peaks = zeros(length(fields), 1);

path2peak_loc_info = fullfile(save_path, 'peak_loc_info');
if ~exist(path2peak_loc_info, 'dir')
    mkdir(path2peak_loc_info);
end

for ii = 1:length(fields)
    field_name = fields{ii};
    simulated_signal = simulation_result.(field_name); simulated_signal = simulated_signal(:);
    real_signal = real_result.(field_name); real_signal = real_signal(:);

    R = corrcoef(simulated_signal, real_signal);
    corr(ii) = R(1,2);
    
    % Compute Cosine similarity
    cosine_similarity = dot(simulated_signal, real_signal) / ...
                        (norm(simulated_signal) * norm(real_signal));
    cos(ii) = cosine_similarity;
    
    % Compute Root Mean Squared Error (RMSE)
    rmse_val = sqrt(mean((simulated_signal - real_signal).^2));
    rmse(ii) = rmse_val;

    if nargin > 2
        % RMSE, but only for the peaks
        fps = 50; cutting_time = 50; peak_prominence = 0.2;
        % pre-compute and store the peak location information
        peak_loc_info_file = fullfile(path2peak_loc_info, strcat(field_name, '.mat'));
        if ~isfile(peak_loc_info_file) % has not been precomputed yet
            [rmse_peaks_val, start_end_indices] = RMSE_for_peaks(real_signal, simulated_signal, peak_prominence, fps, cutting_time);
            save(peak_loc_info_file, "start_end_indices");
        else % use the precomputed start end indices to accelerate
            start_end_indices = load(peak_loc_info_file);
            start_end_indices = start_end_indices.start_end_indices;
            sqr_diff_sum = 0;
            n_element = sum(start_end_indices(2, :) - start_end_indices(1, :)) + size(start_end_indices, 2);
            for jj = 1:size(start_end_indices, 2)
                idx_start = start_end_indices(1, jj);
                idx_end = start_end_indices(2, jj);
                sqr_diff_sum = sqr_diff_sum + sum((real_signal(idx_start:idx_end) - simulated_signal(idx_start:idx_end)).^2);
            end
            rmse_peaks_val = sqr_diff_sum / n_element;
        end
        rmse_peaks(ii) = rmse_peaks_val;
    end
end

if nargin == 2
    metric_vector = [mean(corr); mean(cos); mean(rmse)];
else
    metric_vector = [mean(corr); mean(cos); mean(rmse_peaks)];
end
end