%% do the train-test split, given the struct whose keys are mice names and values are cell arrays of recording path
function cv_splits = createStratifiedTrainTestSplits(datasetStruct, nFolds, rand_seed)
% createStratifiedTrainTestSplits performs stratified k-fold splitting.
%
% Inputs:
%   datasetStruct: struct with mouse names as keys, each containing a cell array of dataset paths.
%   nFolds: number of cross-validation folds (e.g., 5).
%
% Output:
%   cv_splits: struct array (1 x nFolds) with fields:
%       - train: cell array of training dataset paths
%       - test: cell array of testing dataset paths

mice = fieldnames(datasetStruct);

% Initialize the structure for storing splits
cv_splits = struct('train', cell(1, nFolds), 'test', cell(1, nFolds));

% For reproducibility. If you want random, change to rng('shuffle');
rng(rand_seed); 

% First, create stratified splits per mouse
mouseSplits = struct();
for i = 1:length(mice)
    mouse = mice{i};
    datasets = datasetStruct.(mouse);
    n = numel(datasets);
    
    % Shuffle datasets for each mouse
    shuffledIndices = randperm(n);
    datasetsShuffled = datasets(shuffledIndices);
    
    % Divide into approximately equal folds
    foldIndices = crossvalind('Kfold', n, nFolds);
    
    % Store indices per fold
    mouseSplits.(mouse).datasets = datasetsShuffled;
    mouseSplits.(mouse).foldIndices = foldIndices;
end

% Construct the final splits by combining across mice
for fold = 1:nFolds
    trainSet = {};
    testSet = {};
    
    for i = 1:length(mice)
        mouse = mice{i};
        datasetsShuffled = mouseSplits.(mouse).datasets;
        foldIndices = mouseSplits.(mouse).foldIndices;
        
        testSetMouse = datasetsShuffled(foldIndices == fold);
        trainSetMouse = datasetsShuffled(foldIndices ~= fold);
        
        testSet = [testSet; testSetMouse(:)];
        trainSet = [trainSet; trainSetMouse(:)];
    end
    
    cv_splits(fold).train = trainSet;
    cv_splits(fold).test = testSet;
end

end
