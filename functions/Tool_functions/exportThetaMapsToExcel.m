% write theta values to excel
function exportThetaMapsToExcel(thetaMap1, thetaMap2, filename, for_important_fields)
% exportThetaMapsToExcel(thetaMap1, thetaMap2, filename)
%
% Writes the contents of two containers.Map objects to an Excel file.
%   - First column: map keys
%   - Second column: values from thetaMap1
%   - Third column: values from thetaMap2
% The first row will be: {'Key', nameof(thetaMap1), nameof(thetaMap2)}.
%
% Inputs:
%   thetaMap1  containers.Map with scalar or string values
%   thetaMap2  containers.Map with same keys as thetaMap1
%   filename   (optional) string, e.g. 'output.xlsx' (default: 'thetaValues.xlsx')
%
% Example:
%   theta1 = containers.Map({'a','b','c'},{1,2,3});
%   theta2 = containers.Map({'a','b','c'},{10,20,30});
%   exportThetaMapsToExcel(theta1, theta2, 'myThetas.xlsx');

    % Input validation
    assert(isa(thetaMap1,'containers.Map') && isa(thetaMap2,'containers.Map'), ...
        'Both inputs must be containers.Map objects.');
    keys1 = thetaMap1.keys;
    keys2 = thetaMap2.keys;
    assert(isequal(sort(keys1), sort(keys2)), ...
        'thetaMap1 and thetaMap2 must have exactly the same set of keys.');
    
    % Determine output filename
    if nargin < 3 || isempty(filename)
        filename = 'thetaValues.xlsx';
    end

    % Build header row using the variable names
    name1 = inputname(1);  if isempty(name1), name1 = 'theta_manual'; end
    if for_important_fields
        name2 = inputname(2);  if isempty(name2), name2 = 'theta_GA_imp'; end
    else
        name2 = inputname(2);  if isempty(name2), name2 = 'theta_GA_all'; end
    end
    header = {'Key', name1, name2};

    % Sort the keys for consistent ordering
    sortedKeys = sort(keys1);

    % Prepare cell array: (numKeys+1) × 3
    n = numel(sortedKeys);
    outCell = cell(n+1, 3);
    outCell(1, :) = header;

    for i = 1:n
        k = sortedKeys{i};
        outCell{i+1, 1} = k;
        outCell{i+1, 2} = thetaMap1(k);
        outCell{i+1, 3} = thetaMap2(k);
    end

    % Write to Excel (uses writecell, R2019a+)
    writecell(outCell, filename, 'Sheet', 1, 'Range', 'A1');

    fprintf('Wrote %d keys to "%s" with headers [%s, %s].\n', ...
        n, filename, name1, name2);
end