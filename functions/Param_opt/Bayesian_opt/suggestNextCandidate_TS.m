function candMap = suggestNextCandidate_TS(population, cachedFitness, ...
                                           direction_indicator,    ...
                                           param_ranges, refPoint, selected_fields, ...
                                           numDraw, numCand)
% SUGGESTNEXTCANDIDATE_TS  –  one Thompson‑sampling BO step
%
% population         : 1×n cell   – containers.Map of hyper‑parameters
% cachedFitness      : 1×n cell   – m×1 vectors of objective values
% direction_indicator: m×1 (1=minimise, -1=maximise)
% param_ranges       : struct     – each field = [lo hi] row‑vector
% refPoint           : m×1 vector – HV reference (dominated by all)
% selected_fields    : 1xd cell   - field names to be optimized
% numDraw            : scalar     – # TS function draws  (default 1)
% numCand            : scalar     – # Sobol candidate pts (default 1024)
%
% returns candMap    : containers.Map with suggested setting

if nargin < 7, numDraw = 1;    end
if nargin < 8, numCand = 1024; end

% ---------- 1. collect data matrices ------------------------------------
[X, paramNames] = pop2mat(population, selected_fields);          % X: n×d; paramNames: 1xd
Y      = cell2mat(cachedFitness);               % m×n  (each column a point)
Y      = real(Y);
Yadj   = bsxfun(@times, Y, direction_indicator);% all‑minimisation view
m      = size(Y,1);     %# objectives
d      = size(X,2);     %# hyper‑parameters

% ---------- 2. fit one GP per objective ---------------------------------
gp = cell(1,m);
for j = 1:m
    gp{j} = fitrgp(X, Yadj(j,:).', ...
        'BasisFunction','constant', ...
        'KernelFunction','ardsquaredexponential', ...
        'Standardize',true);
end

% ---------- 3. Sobol candidate set --------------------------------------
sob   = sobolset(d,'Skip',1e3,'Leap',1e2);
Xcand = net(sob, numCand);                      % in [0,1]^d

% --- NEW : turn unit cube into actual ranges using param_ranges struct ---
lb = zeros(1,d);  ub = zeros(1,d);
for k = 1:d
    range    = param_ranges.(paramNames{k});    % [lo hi]
    lb(k)    = range(1);
    ub(k)    = range(2);
end
scale  = ub - lb;
Xcand  = bsxfun(@plus, lb, bsxfun(@times, Xcand, scale));  % numCand×d

% ---------- 4. Thompson‑sample & HVI loop -------------------------------
bestGain   = -Inf;
bestVector = [];
hvBase     = hypervolumeND(nondominated(Yadj), refPoint);

for draw = 1:numDraw
    Ydraw = zeros(numCand, m);
    for j = 1:m
        [mu, s2]   = predict(gp{j}, Xcand);
        Ydraw(:,j) = mu + sqrt(max(s2,0)).*randn(numCand,1);
    end
    YdrawAdj = bsxfun(@times, Ydraw, direction_indicator.'); % minimise

    for i = 1:numCand
        newPareto = nondominated([Yadj, YdrawAdj(i,:).']);
        hvNew     = hypervolumeND(newPareto, refPoint);
        gain      = hvNew - hvBase;
        if gain > bestGain
            bestGain   = gain;
            bestVector = Xcand(i,:);
        end
    end
end

% ---------- 5. pack winner into containers.Map --------------------------
candMap = containers.Map;
for k = 1:d
    candMap(paramNames{k}) = bestVector(k);
end
end
