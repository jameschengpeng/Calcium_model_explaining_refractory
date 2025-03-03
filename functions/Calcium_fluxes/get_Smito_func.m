%% modelling the calcium flux from cytosol to mitochondria
function Smito_func = get_Smito_func(K_smito)
Smito_func = @(c_cyto, c_mito) get_Smito(c_cyto, c_mito, K_smito);
end

function Smito = get_Smito(c_cyto, c_mito, K_smito)
if c_mito <= 5 && c_cyto >= 0.2
    Smito = Hill_func(c_cyto, 2, K_smito);
else
    Smito = 0;
end
end