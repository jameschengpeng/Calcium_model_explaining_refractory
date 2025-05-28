%% modelling the calcium flux from cytosol to ER
% return a function handle with c_cyto as the input
function Fserca_func = get_Fserca_func(K_serca, other_settings)
Fserca_func = @(c_cyto, c_ER, init_c_ER) get_Fserca(c_cyto,c_ER, K_serca, other_settings, init_c_ER);
end

function Fserca = get_Fserca(c_cyto, c_ER, K_serca, other_settings, init_c_ER)
scaling_coeff = 1.2;
if ~other_settings('IP3R2') && other_settings('glutamate') % IP3R2-/- but have glutamate transporter
    if c_cyto <= 0.05 || c_ER >= init_c_ER * scaling_coeff
        Fserca = 0.1 * Hill_func(c_cyto, 2, K_serca);% 0.1 micro molar of calcium is a typical value
    else
        % for this modelling, see page 108 of the book computational
        % glioscience, but the dissociation constant is ambiguous
        Fserca = Hill_func(c_cyto, 2, K_serca);
    end
else % normal case
    if c_cyto <= 0.05 || c_ER >= init_c_ER * scaling_coeff
        Fserca = 0;
    else
        % for this modelling, see page 108 of the book computational
        % glioscience, but the dissociation constant is ambiguous
        Fserca = Hill_func(c_cyto, 2, K_serca);
    end
end
end