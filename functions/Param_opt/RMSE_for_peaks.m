%% To evaluate RMSE ONLY on peaks
function [rmse_peaks, start_end_indices] = RMSE_for_peaks(real_signal, pred_signal, peak_prominence, fps, cutting_time)
% filter the real signal to mitigate influence of noise
real_signal = imgaussfilt3(real_signal, fps);

% compute the baseline value
[~, baseline] = response_processing(real_signal, fps, cutting_time);

% make sure all in column vectors
real_signal = real_signal(:);
pred_signal = pred_signal(:);
baseline = baseline(:);

% find the peaks
[pks, locs] = findpeaks(real_signal, 'MinPeakProminence', peak_prominence);

% for each peak, extract the time indices that it reaches 10% amplitude
% build peak only signals
real_peak_only = zeros(size(real_signal));
pred_peak_only = zeros(size(pred_signal));
appending_start = 1;
accept_thres = 0.03;

start_end_indices = zeros(2, length(pks)); % store the start and end indices for each peak (cut at 50% amplitude)

for ii = 1:length(pks)
    loc = locs(ii);
    pk = pks(ii);
    amp = pk - baseline(loc);
    base_amp = baseline + 0.5 * amp;
    
    % indices that signal reaches 10% amplitude
    idx_start = findIntersection(real_signal(1:loc), base_amp(1:loc), accept_thres, 'last'); idx_start = max(1, idx_start);
    if isempty(idx_start)
        idx_start = 1;
    end
    idx_end = findIntersection(real_signal(loc:end), base_amp(loc:end), accept_thres, 'first'); idx_end = min(length(real_signal), idx_end + loc);
    if isempty(idx_end)
        idx_end = length(real_signal);
    end
    len = idx_end - idx_start + 1;

    % append to real_peak_only and pred_peak_only
    real_peak_only(appending_start:(appending_start+len-1)) = real_signal(idx_start:idx_end);
    pred_peak_only(appending_start:(appending_start+len-1)) = pred_signal(idx_start:idx_end);
    appending_start = appending_start+len;
    % update start and end indices
    start_end_indices(1, ii) = idx_start;
    start_end_indices(2, ii) = idx_end;
end

% keep the sub-seq up to appending_start
real_peak_only = real_peak_only(1:appending_start);
pred_peak_only = pred_peak_only(1:appending_start);

rmse_peaks = sqrt(mean((real_peak_only - pred_peak_only).^2));
end

%% for two time series sig1 and sig2, find the first or last intersection point. Two points intersect if they differ no more than accept_thres
% sig1 and sig2: column vectors of same length
% accept_thres: decimal number
% mode: 'first' or 'last'
function idx = findIntersection(sig1, sig2, accept_thres, mode)
% findIntersection  First or last index where two signals intersect.
%
%   idx = findIntersection(sig1, sig2, accept_thres)
%   returns the first index where abs(sig1 - sig2) <= accept_thres.
%
%   idx = findIntersection(sig1, sig2, accept_thres, mode)
%   where mode is 'first' or 'last', returns that respective index.
%
%   If no index satisfies the threshold, idx is returned empty.
%
%   Inputs:
%     sig1, sig2      Vectors of the same length.
%     accept_thres    Nonnegative scalar: intersection threshold.
%     mode            (optional) 'first' (default) or 'last'.
%
%   Output:
%     idx             Scalar index of the intersection, or [].

    % --- Input checking ---
    if nargin<4, mode = 'first'; end
    validateattributes(sig1, {'numeric'}, {'vector'}, mfilename, 'sig1', 1);
    validateattributes(sig2, {'numeric'}, {'vector'}, mfilename, 'sig2', 2);
    validateattributes(accept_thres, {'numeric'}, {'scalar','>=',0}, mfilename, 'accept_thres', 3);
    validatestring(mode, {'first','last'}, mfilename, 'mode', 4);
    
    if numel(sig1)~=numel(sig2)
        error('sig1 and sig2 must have the same number of elements.');
    end
    
    % --- Compute where they "intersect" within threshold ---
    diff_abs = abs(sig1(:) - sig2(:));      % force column vectors
    inds    = find(diff_abs <= accept_thres);
    
    % --- Select first or last as requested ---
    if isempty(inds)
        idx = [];      % no intersection
    else
        switch mode
            case 'first'
                idx = inds(1);
            case 'last'
                idx = inds(end);
        end
    end
end
