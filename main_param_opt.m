clear;
clc;

base_path = fileparts(which('main_param_opt.m'));
data_path = fullfile(base_path, 'Data');

% PLEASE CREATE AND PASTE YOUR OWN SAVE PATH HERE. The train test split information and optimization results
% will be saved in save_path
save_path = '/work/Jamespeng/Astrocyte/Code_for_publish';  % <-- Set this to your own path, e.g., 'C:\Users\YourName\Documents\GA_Optimization_Results'

% Add all subdirectories to the MATLAB path
addpath(genpath(base_path));

%% Hyper parameters for genetic algorithm. Adjust if needed
pop_size = 10; % the population size, divisible by 2
tournament_size = 5; % tournament size for selection
max_generations = 15;
mutation_rate = 0.1;

%% create the savings for train test split information
rand_seed = 21; % <-- Set this to your preferred seed
datasetStruct = build_datasetStruct(data_path);
cv_splits = createStratifiedTrainTestSplits(datasetStruct, 5, rand_seed);
save_path_seed = fullfile(save_path, strcat('seed_', num2str(rand_seed)));
if ~exist(save_path_seed, 'dir')
    mkdir(save_path_seed);
end
save(fullfile(save_path_seed, "cv_split_info.mat"), "cv_splits");
for ii = 1:numel(cv_splits)
    split_folder = fullfile(save_path_seed, strcat('split_', num2str(ii)));
    if ~exist(split_folder, 'dir')
        mkdir(split_folder);
    end
end

%% pre-settings for GA search
fps_org = 30;
fps_upsampled = 50;
other_settings = other_settings_configuration();

% load the manually-tuned theta as the starting point of optimization
theta_init = load(fullfile(base_path, 'theta_init_for_param_opt.mat'));
theta_init = theta_init.theta_init_for_param_opt;

% the important parameters
all_fields = keys(theta_init);  
param_ranges = struct();
% save the bounds to two lists
L_bound = zeros(1, numel(all_fields));
U_bound = zeros(1, numel(all_fields));

for ff = 1:length(all_fields)
    init_value = theta_init(all_fields{ff});
    if strcmp(all_fields{ff}, 'tau_cpkc_degrade')
        lower_bound = 5;
        upper_bound = 100;
    elseif strcmp(all_fields{ff}, 'c_cyto_peak')
        lower_bound = 0.3;
        upper_bound = 1;
    else
        lower_bound = 0.1 * init_value; % make lower bound as 10% of its init value
        upper_bound = 10 * init_value; % make upper bound as 10x of its init value
    end
    param_ranges.(all_fields{ff}) = [lower_bound, upper_bound];
    L_bound(ff) = lower_bound;
    U_bound(ff) = upper_bound;
end

%% train on one split
split_num = 1; % <-- Adjust this split number from 1-5
save_path_split = fullfile(save_path_seed, strcat('split_', num2str(split_num)));
trainSet = cv_splits(split_num).train; testSet = cv_splits(split_num).test;

%% Pre-compute the behavior information to accelerate future computation
path2behavior = fullfile(save_path, 'compiled_behavioral_info');
if ~exist(path2behavior, 'dir')
    mkdir(path2behavior);
end
behavioral_filename = 'behavioral_info.mat';
if ~isfile(fullfile(path2behavior, behavioral_filename))
    allFiles = [trainSet; testSet];
    all_behavior = struct();
    for ii = 1:numel(allFiles)
        aqua_data_path = allFiles{ii};
        [extractedName, extractedLocation, extractedNumber] = extractInfoFromPath(aqua_data_path);
        fieldname = strcat(extractedName, '_', extractedLocation, '_', extractedNumber);
        [resp_seq, run_seq, stim_onset, reward_onset, TT] = ...
                behavioral_data_preprocess(aqua_data_path, fps_org, fps_upsampled);
        all_behavior.(fieldname) = struct();
        all_behavior.(fieldname).resp_seq = resp_seq;
        all_behavior.(fieldname).run_seq = run_seq;
        all_behavior.(fieldname).stim_onset = stim_onset;
        all_behavior.(fieldname).reward_onset = reward_onset;
        all_behavior.(fieldname).TT = TT;
    end
    save(fullfile(path2behavior, behavioral_filename), "all_behavior", "-v7.3");
end

%% perform the evolutionary strategy to optimize the params
% Base theta: a containers.Map with all parameters (fill in the rest as needed)
direction_indicator = [-1; -1; 1];

base_theta = theta_init;
% Initialize population as a cell array of containers.Map objects
population = cell(1, pop_size);
for ii = 1:(pop_size-1)
    candidate = copy_map(base_theta);  % Make a copy of the base theta map
    for j = 1:length(all_fields)
        field = all_fields{j};
        % Randomly sample a value within the bounds for each important parameter.
        candidate(field) = L_bound(j) + (U_bound(j) - L_bound(j)) * rand();
    end
    population{ii} = candidate;
end
population{pop_size} = copy_map(base_theta);
% Evaluate the initial population and cache their fitness.
% Get the current parallel pool.
gen_timer = tic;
pool = gcp();
numCandidates = pop_size;
% Preallocate an array of futures.
futures = parallel.FevalFuture.empty(numCandidates, 0);
% Submit each candidate evaluation as an asynchronous task.
for ii = 1:numCandidates
    futures(ii) = parfeval(pool, @evaluateCandidate_part, 1, population{ii}, other_settings, trainSet, fps_org, fps_upsampled, save_path);
end
% Gather outputs (the fitness metrics) from each future.
cachedFitness = cell(1, numCandidates);
for ii = 1:numCandidates
    cachedFitness{ii} = fetchOutputs(futures(ii));
end
elapsedTime = toc(gen_timer);
savefile = strcat('Generation_', num2str(0), '.mat');
save(fullfile(save_path_split, savefile), "population", "cachedFitness");
disp(['Generation: 0', ' completed in ', num2str(elapsedTime), ' seconds.']);
n_generations = max_generations;
starting_generation = 0;



% loop for the evolution
for generation = 1:n_generations
    % Start timer for the current generation
    gen_timer = tic;
    
    % Form F where each column is an objective vector for a candidate,
    % retrieve from cachedFitness
    F = cell2mat(cachedFitness);
    
    % Perform nondominated sorting
    [ranks, fronts] = nonDominatedSort(F, direction_indicator);
    
    % Select parents via tournament selection (Pareto-aware)
    parents = tournamentSelection(population, ranks, tournament_size, pop_size);
    
    % Create offspring through crossover and mutation
    offspring = cell(1, pop_size);
    for ii = 1:2:(pop_size-1)
        parent1 = parents{ii};
        parent2 = parents{ii+1};
        child1 = crossover(parent1, parent2, all_fields);
        child2 = crossover(parent2, parent1, all_fields);
        offspring{ii} = mutate(child1, param_ranges, all_fields, mutation_rate);
        offspring{ii+1} = mutate(child2, param_ranges, all_fields, mutation_rate);
    end
    
    % Evaluate only offspring
    % Get the current parallel pool.
    pool = gcp();
    % Number of offspring (equals pop_size)
    numCandidates = pop_size;
    % Preallocate an array of futures.
    futures = parallel.FevalFuture.empty(numCandidates, 0);
    % Submit each candidate evaluation as an asynchronous task.
    for ii = 1:numCandidates
        futures(ii) = parfeval(pool, @evaluateCandidate_part, 1, offspring{ii}, other_settings, trainSet, fps_org, fps_upsampled, save_path);
    end
    % Gather outputs (the fitness metrics) from each future.
    offspring_fitness = cell(1, numCandidates);
    for ii = 1:numCandidates
        offspring_fitness{ii} = fetchOutputs(futures(ii));
    end
    % Combine the metrics into the fitness matrix (each column is an objective vector).
    F_offspring = cell2mat(offspring_fitness);  % m x pop_size

    
    % Merge the current population and offspring.
    merged_population = [population, offspring];
    F_merged = [F, F_offspring];  % Combined fitness values: m x (2*pop_size)
    
    % Perform non-dominated sorting on the merged population.
    [merged_ranks, merged_fronts] = nonDominatedSort(F_merged, direction_indicator);
    
    % Select the top candidates based on Pareto rank to form the next generation.
    % We also extract the corresponding fitness values so that we don't need to re-evaluate.
    [population, newFitness] = selectTopCandidates(merged_population, merged_ranks, F_merged, pop_size);
    % convert newFitness (an m x pop_size matrix) into a cell array
    cachedFitness = cell(1, pop_size);
    for ii = 1:pop_size
        cachedFitness{ii} = newFitness(:, ii);
    end

    % Save the population and finalFitness
    savefile = strcat('Generation_', num2str(generation + starting_generation), '.mat');
    save(fullfile(save_path_split, savefile), "population", "cachedFitness");
    
    % Measure the elapsed time and display it
    elapsedTime = toc(gen_timer);
    disp(['Generation: ', num2str(generation + starting_generation), ' completed in ', num2str(elapsedTime), ' seconds.']);
end

