function HV = hypervolumeND(F, r)
% HYPERVOLUMEND  Exact hypervolume of a non‑dominated set.
%
%   HV = hypervolumeND(F, r)
%
%   F : d×n matrix with objective columns (to be MINIMISED)
%   r : d×1 reference point that is dominated by every point in F
%
%   Removes dominated points automatically, then applies the
%   Fonseca‑Paquete recursive decomposition.
%
%   (c) 2025  —  feel free to use / adapt

    %----- 0. Input hygiene -----------------------------------------------
    if size(F, 1) ~= numel(r)
        error('Row‑dimension of F (%d) must equal length of r (%d).', size(F,1), numel(r));
    end
    F = F(:, all(F < r, 1));       % discard any point outside the ref box
    if isempty(F)
        HV = 0.0;  return;
    end

    %----- 1. Pareto filter (non‑dominated)---------------------------------
    F = nondominated(F);           % custom function below, O(n·d)

    %----- 2. Recursive HV -------------------------------------------------
    HV = hvRecursive(F, r);
end

%% ==========================================================================

function P = nondominated(F)
% Fast non‑dominated filter for d<=10 and n <= 1e4
    [d, n] = size(F);
    keep = true(1,n);
    for i = 1:n-1
        if ~keep(i), continue, end
        fi = F(:, i);
        % point j dominates i  <=>  F(:,j) <= fi & any(F(:,j) < fi)
        dom = all(F(:, i+1:n) <= fi, 1) & any(F(:, i+1:n) < fi, 1);
        keep(i+find(dom)) = false;
        % if i is dominated, mark and break inner loop early
        if any(all(fi >= F(:, i+1:n) ,1) & any(fi > F(:, i+1:n),1))
            keep(i) = false;
        end
    end
    P = F(:, keep);
end

%% ==========================================================================

function hv = hvRecursive(F, r)
% Recursive hypervolume (Fonseca, Paquete & López‑Ibañez 2006)
    [d, n] = size(F);

    % ---- Base cases -----------------------------------------------------
    if n == 0
        hv = 0.0;
        return
    elseif d == 1
        hv = r - min(F);
        return
    elseif d == 2
        hv = hv2D(F, r);           % efficient strip algorithm below
        return
    end

    % ---- d >= 3 : recursive decomposition -------------------------------
    % sort by first objective (ascending => best first)
    [~, idx] = sort(F(1, :), 'ascend');
    F = F(:, idx);

    hv   = 0.0;
    prev = r(1);
    for j = n:-1:1                 % iterate worst → best
        w = prev - F(1, j);        % width in objective‑1
        if w > 0
            slice = F(2:end, 1:j); % points with smaller/equal obj‑1
            hv   = hv + w * hvRecursive(slice, r(2:end));
            prev = F(1, j);
        end
    end
end

%% ==========================================================================

function area = hv2D(F, r)
% Exact 2‑D hypervolume (area of rectangular strips)
% F : 2×n (non‑dominated, minimisation)
    [~, idx] = sort(F(1, :), 'ascend');
    F = F(:, idx);

    area   = 0.0;
    prev_x = r(1);
    for j = size(F,2):-1:1
        width  = prev_x - F(1,j);
        height = r(2)  - min(F(2,1:j));
        area   = area + width * height;
        prev_x = F(1, j);
    end
end
