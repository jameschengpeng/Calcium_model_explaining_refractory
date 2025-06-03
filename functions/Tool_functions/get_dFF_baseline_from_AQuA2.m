%% Main function: get_dFF_baseline_from_AQuA2
function [dFF_list, dFF_list_by_aqua, baseline_list] = get_dFF_baseline_from_AQuA2(AQuA2_filepath)
    % Set initial number of workers and maximum number of retries.
    maxRetries = 5;
    currentWorkers = 26;  % starting with 26 workers

    % Load file on the client to determine number of events.
    S_client = load(AQuA2_filepath);
    fields = fieldnames(S_client);
    matchingFields = fields(startsWith(fields, 'res'));
    res_dbg_client = S_client.(matchingFields{1});
    n_evts = length(res_dbg_client.evt1);

    % Preallocate output arrays.
    dFF_list = zeros(n_evts, 1);
    dFF_list_by_aqua = zeros(n_evts, 1);
    baseline_list = zeros(n_evts, 1);

    for attempt = 1:maxRetries
        try
            % (Re)start a parallel pool with the current number of workers.
            pool = gcp('nocreate');
            if isempty(pool) || pool.NumWorkers ~= currentWorkers
                if ~isempty(pool)
                    delete(pool);
                end
                parpool('local', currentWorkers);
            end

            % Use a parallel constant to load the file once per worker.
            pc = parallel.pool.Constant(@() load(AQuA2_filepath));
            
            % Run the parfor loop.
            parfor evt_num = 1:n_evts
                try
                    % Each worker retrieves its local copy.
                    S_local = pc.Value;
                    localFields = fieldnames(S_local);
                    localMatchingFields = localFields(startsWith(localFields, 'res'));
                    local_res_dbg = S_local.(localMatchingFields{1});

                    % Process the event.
                    [average_baseline, max_dFF, max_dFF_by_aqua] = get_baseline_dFF(local_res_dbg, evt_num);
                    dFF_list(evt_num) = max_dFF;
                    dFF_list_by_aqua(evt_num) = max_dFF_by_aqua;
                    baseline_list(evt_num) = average_baseline;
                catch innerME
                    % Rethrow with a custom identifier so the outer catch can detect it.
                    error('ParallelIterationFailed:WorkerError', ...
                        'Error in parfor iteration %d: %s', evt_num, innerME.message);
                end
            end

            % If the loop completes successfully, exit the retry loop.
            fprintf('Processing completed using %d workers.\n', currentWorkers);
            return;
        catch ME
            % Check for out-of-memory error (case-insensitive check)
            if contains(lower(ME.message), 'out of memory') || contains(lower(ME.message), 'memory allocation')
                fprintf('Out of Memory error encountered with %d workers. Retrying with half the workers...\n', currentWorkers);
                % Halve the number of workers, ensuring at least one remains.
                currentWorkers = max(floor(currentWorkers / 2), 1);
                % Delete any existing pool before retrying.
                pool = gcp('nocreate');
                if ~isempty(pool)
                    delete(pool);
                end
            else
                % For other errors, rethrow the exception.
                rethrow(ME);
            end
        end
    end

    % If we exit the loop without returning, then we failed to process after several attempts.
    error('Failed to complete processing after %d attempts due to memory constraints.', maxRetries);
end

%% Helper function: get_baseline_dFF
function [average_baseline, max_dFF, max_dFF_by_aqua] = get_baseline_dFF(res_dbg, evt_num)
    [evt_dFF, baseline, unique_z] = get_evt_dFF(res_dbg.evt1{evt_num}, res_dbg.datOrg1);
    average_baseline = mean(baseline(unique_z));
    max_dFF = max(evt_dFF(unique_z));
    max_dFF_by_aqua = max(squeeze(res_dbg.dffMat1(evt_num, unique_z, 1)));
end

%% Helper function: uniqueCoordinates
function [XX, YY] = uniqueCoordinates(X, Y)
    % Ensure X and Y are column vectors
    if ~iscolumn(X)
        X = X';
        Y = Y';
    end
    points = [X Y];
    % Find unique coordinate pairs
    uniquePoints = unique(points, "row");
    % Extract unique X and Y coordinates
    XX = uniquePoints(:, 1);
    YY = uniquePoints(:, 2);
end

%% Helper function: get_evt_dFF
function [evt_dFF, baseline, unique_z] = get_evt_dFF(event, data)
    data = squeeze(data); % data is 4D, make it 3D
    [x, y, z] = ind2sub(size(data), event);
    binary = zeros(size(data));
    binary(event) = 1;
    [XX, YY] = uniqueCoordinates(x, y);
    ind2D = sub2ind([size(data,1), size(data,2)], XX, YY);
    rough_mean = zeros(1, size(data, 3));
    for ii = 1:length(rough_mean)
        data_slice = data(:, :, ii);
        rough_mean(ii) = mean(data_slice(ind2D));
    end
    [~, baseline] = response_processing(rough_mean, 3, 200);
    unique_z = unique(z);
    event_mean = zeros(1, size(data, 3)); 
    for ii = 1:length(unique_z)
        zz = unique_z(ii);
        binary_slice = binary(:, :, zz);
        data_slice = data(:, :, zz);
        event_mean(zz) = mean(data_slice(binary_slice == 1));
    end
    dF = max(baseline, event_mean) - baseline;
    evt_dFF = dF ./ baseline;
end
