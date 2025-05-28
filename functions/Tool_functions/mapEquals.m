function tf = mapEquals(m1, m2)
%MAPEQUALS True if two containers.Map objects are equal in every attribute.
%
%   TF = MAPEQUALS(M1, M2) returns logical true if
%   • M1 and M2 are both containers.Map objects
%   • They have the same KeyType, ValueType, and Count
%   • They contain the exact same set of keys (order-independent)
%   • Corresponding values are equal, using ISEQUALN
%
%   Nested containers.Map values are handled recursively.

    % 1) Basic class/type checks
    if ~isa(m1, 'containers.Map') || ~isa(m2, 'containers.Map')
        tf = false;
        return
    end

    % 2) Simple property comparisons
    if m1.Count      ~= m2.Count  || ...
       ~strcmp(m1.KeyType,   m2.KeyType) || ...
       ~strcmp(m1.ValueType, m2.ValueType)
        tf = false;
        return
    end

    % 3) Compare key sets (order-independent)
    k1 = sort(m1.keys);
    k2 = sort(m2.keys);
    if ~isequal(k1, k2)
        tf = false;
        return
    end

    % 4) Compare each value
    for i = 1:numel(k1)
        v1 = m1(k1{i});
        v2 = m2(k1{i});
        if isa(v1, 'containers.Map') && isa(v2, 'containers.Map')
            % Recursive comparison for nested maps
            if ~mapEquals(v1, v2)
                tf = false;
                return
            end
        else
            if ~isequaln(v1, v2)
                tf = false;
                return
            end
        end
    end

    % All tests passed
    tf = true;
end
