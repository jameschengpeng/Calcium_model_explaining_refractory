function [newPath, isNative] = remapPath(oldPath, dstRoot, anchor)
% REMAPPATH  Cross-platform path fixer with automatic OS detection.
%
%   [newPath, isNative] = REMAPPATH(oldPath, dstRoot, anchor)
%
%   oldPath   : char | string | cell array   – absolute path(s) from any OS
%   dstRoot   : char | string  (optional)   – where to re-root *if* conversion
%               is needed.  Example: 'C:\Users\james\CBIL'  or '/work/Jamespeng'
%   anchor    : char | string  (optional)   – folder that marks the start of
%               the relative part to keep.  Example: 'Astrocyte'
%
%   newPath   : same shape as oldPath, with separators native to *this* MATLAB
%   isNative  : logical(s) indicating whether each oldPath already matched
%               the host OS (no remap needed)
%
%   ─────────────────────────────────────────────────────────────────────────
%   • If dstRoot & anchor are omitted, the function just normalises the
%     separators and tells you whether the original was native.
%   • Works recursively on cell arrays.
%   • For mismatched paths **and** provided dstRoot+anchor, the relative
%     sub-path (anchor → end) is glued onto dstRoot via fullfile.
%
%   Examples
%   --------
%   1) Just check/clean a single path
%       p = '/home/james/data/file.mat';
%       [fixed,native] = remapPath(p);
%
%   2) Convert Linux → Windows
%       src = '/work/Jamespeng/Astrocyte/Exp1/file.mat';
%       dst = 'C:\Users\james\CBIL';
%       new = remapPath(src,dst,'Astrocyte');
%
%   3) Batch-convert Windows → Linux
%       list = {'C:\Users\james\CBIL\Astrocyte\Rec1\a.mat'
%               'C:\Users\james\CBIL\Astrocyte\Rec2\b.mat'};
%       linuxRoot = '/work/Jamespeng';
%       [newList,~] = remapPath(list,linuxRoot,'Astrocyte');

    % ── Batch recursion ───────────────────────────────────────────────────
    if iscell(oldPath)
        [newPath, isNative] = cellfun(@(p) remapPath(p,dstRoot,anchor), ...
                                      oldPath, 'UniformOutput', false);
        isNative = cell2mat(isNative);            % flatten logical array
        return
    end

    % ── 1. What style is the *source* path? ───────────────────────────────
    srcStyle = detectStyle(oldPath);              % 'windows', 'posix', 'unknown'

    % ── 2. What style does *this* MATLAB expect? ──────────────────────────
    if ispc
        hostStyle = 'windows';
    else
        hostStyle = 'posix';                      % covers Linux & macOS
    end

    isNative = strcmp(srcStyle,hostStyle);

    % ── 3. If native, just pretty-print the separators and return ─────────
    if isNative || nargin < 3                     % (no anchor ⇒ no remap)
        newPath = fixSeps(oldPath, filesep);
        return
    end

    % ── 4. Otherwise we need to re-root the path ──────────────────────────
    assert(~isempty(dstRoot) && ~isempty(anchor), ...
          ['Path style mismatch (' srcStyle ' → ' hostStyle ') and no ', ...
           'dstRoot+anchor supplied.']);

    newPath = actuallyRemap(oldPath,dstRoot,anchor);
end

% ======================================================================== %
function style = detectStyle(p)
    p = char(p);                                  % be sure it's 1×n char
    if ~isempty(regexp(p,'^[A-Za-z]:[\\/]', 'once')) || startsWith(p,'\\')
        style = 'windows';
    elseif startsWith(p,'/')
        style = 'posix';
    else
        style = 'unknown';
    end
end

function out = fixSeps(p,sep)
    p = strrep(p,'\','/');                        % make uniform
    if sep == '\', p = strrep(p,'/','\'); end     % convert if needed
    out = p;
end

function newPath = actuallyRemap(oldPath,dstRoot,anchor)
    % --- unify separators so string ops are predictable
    src = strrep(char(oldPath),'\','/');
    dst = strrep(char(dstRoot),'\','/');

    needle = ['/' anchor '/'];
    idx    = strfind(src,needle);
    assert(~isempty(idx), ['Anchor "' anchor '" not found in "' oldPath '".']);

    rel    = src(idx+1:end);                      % drop leading '/'
    parts  = strsplit(rel,'/');

    newPath = fullfile(dst,parts{:});             % host-native sep
end
