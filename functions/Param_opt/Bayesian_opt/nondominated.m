function [Fnd, keep] = nondominated(F)
% NONDOMINATED   Return the non‑dominated subset of a point matrix.
%
%   [Fnd, keep] = nondominated(F)
%
%   F     : d×n numeric matrix (each column is an objective vector,
%           and **smaller is better** for every objective).
%
%   Fnd   : d×k matrix containing only the non‑dominated columns (k ≤ n)
%   keep  : 1×n logical mask, true for non‑dominated points
%
%   A point p dominates q  ⇔  p ≤ q component‑wise  AND  p < q in at least one
%   component.  Complexity 𝒪(n²·d) — fine for n up to a few thousand.

    [d, n] = size(F);
    keep   = true(1, n);

    for i = 1:n
        if ~keep(i), continue, end

        % Any point that is element‑wise ≤ the i‑th *and* strictly < in one
        % component dominates the i‑th.  Mark such points as dominated.
        dominatedByI = all(bsxfun(@le, F(:,i), F), 1) & ...
                       any(bsxfun(@lt, F(:,i), F), 1);
        keep(dominatedByI) = false;

        % If the current point i is itself dominated by **any** other point,
        % mark it and skip the expensive comparisons later.
        if any(all(bsxfun(@le, F, F(:,i)),1) & ...
               any(bsxfun(@lt, F, F(:,i)),1))
            keep(i) = false;
        end
    end

    Fnd = F(:, keep);
end
