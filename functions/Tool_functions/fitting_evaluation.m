%% perform multi metric evaluation
function [corr, cosine_similarity, rmse] = fitting_evaluation(filepath, resp_seq)
result = load(filepath);
chemical_table = result.chemical_table;
c_cyto = chemical_table{:,"cytosolic_Ca"}; simulated_signal = c_cyto(:);
real_signal = resp_seq(:);

simulated_signal = scale_resp(simulated_signal, 0.05, 1);
real_signal = scale_resp(real_signal, 0.05, 1);

R = corrcoef(simulated_signal, real_signal);
corr = R(1,2);

cosine_similarity = dot(simulated_signal, real_signal) / ...
                    (norm(simulated_signal) * norm(real_signal));

rmse = sqrt(mean((simulated_signal - real_signal).^2));
end