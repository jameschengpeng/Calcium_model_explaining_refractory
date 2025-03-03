%% modelling the calcium flux from ER to cytosol
% c_cyto and ip3 affects the probability of binding, i.e., IP3R binds with
% both calcium and ip3
% atp affects the probability of opening given that IP3R is binded with
% calcium and ip3
% return a function handle for the calcium flux from ER to cytosol, with
% c_cyto, c_ER, ATP as the variables, and IP3 & n (num of IP3R2) as constants
function [FER_func, c_cyto_effect, ATP_effect, IP3_effect] = get_FER_func(theta, K_concentration, ROS_cyto, act_cPKC_particle, other_settings)
IP3R1_percentage = theta('IP3R1_percentage');
IP3R2_percentage = 1 - IP3R1_percentage;
%% parameters for IP3R2
mu1 = 0.4195;
sigma1 = 0.6607;
mu2 = 2.1490;
b1 = 3.0918;
sigma2 = 2.1530;
c1 = 1.3039;

%% parameters for IP3R1
b1_R1 = 0.4395;
mu1_R1 = 0.0264;
sigma1_R1 = 0.8146;
c1_R1 = 0.8194;
c2_R1 = 0.5942;

%%
cPKC = mean(act_cPKC_particle(:,2));
% IP3R2
c_cyto_effect = @(c_cyto) exp(-(log10(c_cyto) + mu1)^2/(2*sigma1^2));
ATP_effect = @(ATP, IP3) exp(-(log10(ATP) - mu2 + b1*log10(IP3))^2/(2*sigma2^2));
IP3_effect = @(IP3, ATP) IP3^2/(IP3^2 + (c1-log10(ATP))^2);
% IP3R1
c_cyto_ATP_effect_R1 = @(c_cyto, ATP) exp(-(log10(c_cyto)-(-b1_R1*log10(ATP)+mu1_R1))^2 / (2 * sigma1_R1^2));
IP3_effect_R1 = @(IP3) (log10(IP3)+c2_R1)^2 / ((log10(IP3)+c2_R1)^2 + c1_R1^2);

ROS = ROS_cyto;
ROS_basis_effect = 0.7;
ROS_effect = ROS_basis_effect + (1-ROS_basis_effect) * Hill_func(ROS, 1, 0.3);
PKA_effect = @(PKA_particle) 0.4 * Hill_func(mean(PKA_particle(:,2)), 1, 3);
% n is removed. Useless
FER_func = @(c_cyto, c_ER, ATP, IP3, PKA_particle) ...
    (IP3R2_percentage * (c_cyto_effect(c_cyto) * IP3_effect(IP3, ATP) * ATP_effect(ATP, IP3) * ROS_effect + PKA_effect(PKA_particle)) ...
    + IP3R1_percentage * (c_cyto_ATP_effect_R1(c_cyto, ATP) * IP3_effect_R1(IP3) * ROS_effect + PKA_effect(PKA_particle))) * ...
    concentration_effect(c_cyto, c_ER, K_concentration);
if ~other_settings('IP3R2')
    FER_func = @(c_cyto, c_ER, ATP, IP3) 0;
end
end


% since transfering calcium from ER to cytosol does not consume ATP
% it should depend on the concentration difference in the two compartments
function effect = concentration_effect(c_cyto, c_ER, K_concentration)
if c_cyto >= c_ER
    effect = 0;
else
    diff = c_ER - c_cyto;
    effect = Hill_func(diff, 2, K_concentration);
    % effect = 1;
end
end
