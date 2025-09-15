%% call this function when you wish to modify the theta_configuration.m, using a customized theta (containers.Map)
function syncThetaInitFile(theta, inFile, outFile)
%SYNCTHETAINITFILE Sync theta(...) assignments in inFile to the values in
% the provided containers.Map "theta". Preserves trailing comments.
%
% - Works with function files (inserts appended keys before the final "end")
% - Skips lines whose trailing comment mentions "NOT USED" or "PLEASE SKIP"
% - Appends keys in theta that are missing from the file, with a timestamp

    if nargin < 3, outFile = inFile; end
    if ~isa(theta, 'containers.Map')
        error('First argument must be a containers.Map (your theta).');
    end

    txt   = fileread(inFile);
    lines = regexp(txt, '\r\n|\n|\r', 'split');  % OS-agnostic
    updatedLines = lines;
    seenKeys = containers.Map('KeyType','char','ValueType','logical');

    % Token-based regex (portable across MATLAB versions):
    % Groups:
    %  1: pre  ("   theta('key') = ")
    %  2: key  (inside quotes)
    %  3: val  (old value literal)
    %  4: semi (semicolon/whitespace)
    %  5: comment (optional trailing %... )
    pat = '^(\\s*theta\\(\\s*[''"]([^''"]+)[''"]\\s*\\)\\s*=\\s*)(.*?)(;?\\s*)(%.*)?$';

    for i = 1:numel(lines)
        t = regexp(lines{i}, pat, 'tokens', 'once');
        if isempty(t), continue; end

        pre  = t{1};
        key  = t{2};
        val  = t{3};
        semi = t{4};
        if numel(t) >= 5 && ~isempty(t{5}), comment = t{5}; else, comment = ''; end

        % Skip lines marked NOT USED / PLEASE SKIP
        if ~isempty(comment) && ~isempty(regexpi(comment, 'NOT\s+USED|PLEASE\s+SKIP'))
            seenKeys(key) = true;
            continue; % leave line unchanged
        end

        if isKey(theta, key)
            try
                lit = value2literal(theta(key));
            catch ME
                warning('Key "%s": %s (leaving line unchanged)', key, ME.message);
                lit = val; % fallback: keep original
            end
            if ~isempty(comment)
                updatedLines{i} = [pre, lit, semi, comment];
            else
                updatedLines{i} = [pre, lit, semi];
            end
            seenKeys(key) = true;
        end
    end

    % Append any theta keys not present in file
    toAppend = {};
    allKeys = theta.keys;
    tstamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    for j = 1:numel(allKeys)
        k = allKeys{j};
        if ~isKey(seenKeys, k)
            lit = value2literal(theta(k));
            toAppend{end+1} = sprintf('theta(%s) = %s; %% added by sync on %s', ...
                quoteStr(k), lit, tstamp); %#ok<AGROW>
        end
    end

    % If file has a terminal "end", insert before it
    if ~isempty(toAppend)
        endIdx = find(cellfun(@(L) ~isempty(regexp(L, '^\s*end\s*(%.*)?$', 'once')), updatedLines), 1, 'last');
        if ~isempty(endIdx)
            updatedLines = [updatedLines(1:endIdx-1), toAppend, updatedLines(endIdx:end)];
        else
            updatedLines = [updatedLines, toAppend];
        end
    end

    % Write out
    outTxt = strjoin(updatedLines, newline);
    fid = fopen(outFile, 'w');
    if fid < 0, error('Could not open "%s" for writing.', outFile); end
    fwrite(fid, outTxt);
    fclose(fid);

    fprintf('Synced %s → %s\n', inFile, outFile);
end

function lit = value2literal(v)
    if isnumeric(v)
        lit = mat2str(v, 15);
    elseif islogical(v)
        if isscalar(v), lit = char(string(v)); else, lit = mat2str(v); end
    elseif ischar(v)
        lit = quoteStr(v);
    elseif isa(v,"string")
        if isscalar(v)
            lit = quoteStr(char(v));
        else
            c = cellstr(v);
            inner = cellfun(@quoteStr, c, 'UniformOutput', false);
            lit = ['{', strjoin(inner, ', '), '}'];
        end
    elseif iscell(v)
        if isvector(v)
            parts = cellfun(@value2literal, v, 'UniformOutput', false);
            lit = ['{', strjoin(parts, ', '), '}'];
        else
            parts = cellfun(@value2literal, v(:).', 'UniformOutput', false);
            lit = ['{', strjoin(parts, ', '), '}'];
        end
    else
        error('Unsupported value type: %s', class(v));
    end
end

function q = quoteStr(s)
    if ~ischar(s), s = char(s); end
    q = ['''', strrep(s, '''', ''''''), ''''];
end
