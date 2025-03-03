%% modelling for cAMP
% activated by AC, constant degradation and degradation by
% phosphodiesterase (PDE)
% normally, ATP concentration (in mM) while cAMP concentration (micro M)
% we can assume ATP is infinite when computing cAMP production

% Ni, Qiang, et al. "Signaling diversity of PKA achieved via a Ca2+-cAMP-PKA oscillatory circuit." Nature chemical biology 7.1 (2011): 34-40.

function dcAMP_dt_func = get_dcAMP_dt_func(theta)
dcAMP_dt_func = @(cAMP, AC_particle, PDE_particle) AC_activation(AC_particle, theta) - PDE_deactivation(PDE_particle, cAMP, theta);
end

% normally, ATP concentration (in mM) while cAMP concentration (nM)
% we can assume ATP is infinite when computing cAMP production
% so we do not incorporate Michealis Menten kinetics here
function activate_rate = AC_activation(AC_particle, theta)
b_cAMP_produce = theta('b_cAMP_produce');
AC = mean(AC_particle(:,2));
activate_rate = b_cAMP_produce * AC;
end

% Michealis Menten kinetics
function deactivate_rate = PDE_deactivation(PDE_particle, cAMP, theta)
b_cAMP_degrade1 = theta('b_cAMP_degrade1');
b_cAMP_degrade2 = theta('b_cAMP_degrade2');
PDE = mean(PDE_particle(:,2));
deactivate_rate = b_cAMP_degrade1 * PDE * Hill_func(cAMP, 1, 4) + b_cAMP_degrade2 * cAMP;
end

