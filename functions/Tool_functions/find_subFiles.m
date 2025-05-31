function subFiles = find_subFiles(path)
files = dir(path);
subFiles = files(~[files.isdir]);
end