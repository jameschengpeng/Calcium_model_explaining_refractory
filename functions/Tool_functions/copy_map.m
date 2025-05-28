%% copy a containers.Map, pass by value
function map2 = copy_map(map1)
orgKeys = map1.keys;
orgValues = map1.values;
map2 = containers.Map(orgKeys, orgValues);
end