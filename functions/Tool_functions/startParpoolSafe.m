function p = startParpoolSafe(varargin)
%STARTPARPOOLSAFE Start or reuse a local parpool with conservative defaults.
%
% Usage:
%   p = startParpoolSafe();              % auto workers (cores-1)
%   p = startParpoolSafe(12);            % request 12 workers
%   p = startParpoolSafe('local', 12);   % specify profile and workers
%
% Notes:
%   - Reuses an existing pool when possible.
%   - Uses conservative defaults to avoid worker crashes due to RAM limits.
%   - Falls back to fewer workers if startup fails.

    % ---- Input parsing ----
    profile = 'local';
    reqN = [];

    if nargin == 1
        reqN = varargin{1};
    elseif nargin == 2
        profile = varargin{1};
        reqN = varargin{2};
    elseif nargin > 2
        error('Usage: startParpoolSafe(), startParpoolSafe(N), or startParpoolSafe(profile, N).');
    end

    c = parcluster(profile);
    maxAllowed = c.NumWorkers;

    % ---- Worker selection ----
    if isempty(reqN)
        nCores = feature('numcores');
        numWorkers = max(1, min(maxAllowed, nCores - 1));
    else
        validateattributes(reqN, {'numeric'}, {'scalar','integer','positive'});
        numWorkers = min(reqN, maxAllowed);
    end

    % ---- Reuse existing pool if possible ----
    p0 = gcp('nocreate');
    if ~isempty(p0)
        sameProfile = strcmpi(p0.Cluster.Profile, c.Profile);
        sameSize = (p0.NumWorkers == numWorkers);

        if sameProfile && (isempty(reqN) || sameSize)
            p = p0;
            return;
        end

        delete(p0);
    end

    % ---- Start pool with fallback ----
    try
        p = parpool(c, numWorkers);
    catch ME
        fallback = max(1, min(4, maxAllowed));
        warning(['Failed to start parpool with %d workers (%s).\n' ...
                 'Retrying with %d workers.\nReason: %s'], ...
                 numWorkers, c.Profile, fallback, ME.message);
        p = parpool(c, fallback);
    end
end
