%% Hill function
function h = Hill_func(x, N, c)
h = (x.^N)./(x.^N + c^N);
end