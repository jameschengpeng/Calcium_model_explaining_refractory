%% nonlinear transform the fluorescent signal
function transformed_signal = fluorescent_signal_transform(fluo_signal, n, Kd, min_val, max_val)
    % This function transforms a fluorescence signal into a calcium concentration
    % time series and scales it to a specified range [min_val, max_val].
    %
    % Parameters:
    %   fluo_signal - Vector of fluorescence intensity values
    %   n           - Hill coefficient
    %   Kd          - Dissociation constant (in the same units as the desired output)
    %   min_val     - Desired minimum value of the transformed signal
    %   max_val     - Desired maximum value of the transformed signal
    %
    % Returns:
    %   transformed_signal - Vector of transformed and scaled calcium concentrations

    % Ensure the fluorescence signal is a column vector
    fluo_signal = fluo_signal(:);

    % Determine the minimum and maximum fluorescence intensities
    F_min = min(fluo_signal);
    F_max = max(fluo_signal);

    % Check for potential division by zero
    if F_max == F_min
        error('Fluorescence signal has no variation; cannot perform transformation.');
    end

    % Transform fluorescence to calcium concentration using the inverse Hill equation
    Ca_concentration = Kd * ((fluo_signal - F_min) ./ (F_max - fluo_signal + 0.01)).^(1 / n);

    % Determine the minimum and maximum calcium concentrations
    Ca_min = min(Ca_concentration);
    Ca_max = max(Ca_concentration);

    % Check for potential division by zero
    if Ca_max == Ca_min
        error('Calcium concentration has no variation; cannot perform scaling.');
    end

    % Scale the calcium concentration to the desired range [min_val, max_val]
    transformed_signal = min_val + ((Ca_concentration - Ca_min) / (Ca_max - Ca_min)) * (max_val - min_val);
end
