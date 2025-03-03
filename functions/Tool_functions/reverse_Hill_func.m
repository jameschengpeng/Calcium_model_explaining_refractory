%% reverse of Hill function
function h = reverse_Hill_func(x, N, c)
h = c^N/(x^N + c^N);
end