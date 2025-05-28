function refPoint = buildRefPoint(cachedFitness, direction_indicator, margin)
% BUILDREFPOINT   Choose a reference point for hyper‑volume.
%
%   refPoint = buildRefPoint(cachedFitness, direction_indicator, margin)
%
%   cachedFitness      : 1×n cell, each cell m×1 column with the metrics
%   direction_indicator: m×1   (1 = minimise,  -1 = maximise)
%   margin             : scalar (0.10 ⇒ 10 % worse than worst) – optional
%
%   refPoint           : m×1 column in **minimisation** coordinates

    if nargin < 3,  margin = 0.10;  end

    Y     = cell2mat(cachedFitness);                 % m×n
    Yadj  = bsxfun(@times, Y, direction_indicator);  % cast → minimise
    m     = size(Yadj,1);

    refPoint = zeros(m,1);
    for j = 1:m
        worst = max(Yadj(j,:));      % larger = worse (in minimisation view)
        if worst >= 0
            % e.g. RMSE in [0,1] or any loss‑like metric
            refPoint(j) = worst * (1 + margin);
        else
            % e.g. −(correlation) or −(cosineSim), values are negative
            refPoint(j) = worst + margin * abs(worst);
        end
    end
end
