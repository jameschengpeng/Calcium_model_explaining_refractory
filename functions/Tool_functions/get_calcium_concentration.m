%% Estimating the calcium concentration during calcium events
% take vector inputs
function calcium_concentration = get_calcium_concentration(dFF, baseline, n, Kd, max_min_ratio, ca_base)
Fmin = get_Fmin(n, Kd, max_min_ratio, ca_base, baseline);
calcium_concentration = get_adjusted_calcium(n, Kd, max_min_ratio, dFF, baseline, Fmin);
end

function Fmin = get_Fmin(n, Kd, max_min_ratio, ca_base, F0)
hill_value = Hill_func(ca_base, n, Kd);
Fmin = F0 ./ (1 + hill_value*(max_min_ratio-1));
end

function adjusted_calcium = get_adjusted_calcium(n, Kd, max_min_ratio, dFF, F0, Fmin)
hill_value = ((dFF + 1).*F0 - Fmin) ./ ((max_min_ratio-1).*Fmin); % this is a vector
hill_value = min(hill_value, 0.99); % however, sometimes it can be > 1, which is because max_min_ratio is at lower bound
ca_power_n = (Kd^n .* hill_value) ./ (1 - hill_value);
adjusted_calcium = ca_power_n .^ (1/n);
end