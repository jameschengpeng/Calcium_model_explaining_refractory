%% find subdirectories
function subDirs = find_subDirs(path)
allFiles = dir(path);
isSubDir = [allFiles.isdir];
subDirs = {allFiles(isSubDir).name};
subDirs = subDirs(~ismember(subDirs,{'.','..'}));
end