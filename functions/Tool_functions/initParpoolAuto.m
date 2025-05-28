function initParpoolAuto()

    % Close existing pool if active
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end

    % Default fallback
    numWorkers = 26;
    c = parcluster('local');
    maxAllowed = c.NumWorkers;

    if isunix
        % Assume Linux server
        numWorkers = 96;

    elseif ispc
        % Get CPU name via PowerShell
        try
            [~, cpuInfo] = system('powershell -Command "Get-CimInstance Win32_Processor | Select-Object -ExpandProperty Name"');
            cpuInfo = strtrim(cpuInfo);

            if contains(cpuInfo, "Xeon", 'IgnoreCase', true)
                numWorkers = 26;
            elseif contains(cpuInfo, "i7-14700F", 'IgnoreCase', true)
                numWorkers = 16;
            else
                warning('Unrecognized Windows CPU: "%s". Using default of 4 workers.', cpuInfo);
            end
        catch
            warning('Failed to retrieve CPU info using PowerShell. Using default of 4 workers.');
        end
    end

    % Enforce max limit
    if numWorkers > maxAllowed
        warning('Requested %d workers exceeds cluster limit of %d. Using max allowed.', numWorkers, maxAllowed);
        numWorkers = maxAllowed;
    end

    % Start the pool
    parpool('local', numWorkers);
    fprintf('Started parpool with %d workers on %s\n', numWorkers, getHostName());
end

function hostname = getHostName()
    if ispc
        [~, hostname] = system('hostname');
    else
        [~, hostname] = system('uname -n');
    end
    hostname = strtrim(hostname);
end
