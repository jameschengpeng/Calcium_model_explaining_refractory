%% Function to plot vertical lines and create a handle for the legend
function h = plotVerticalLines(x, indices, color)
    if isempty(indices)
        h = []; % No handle if the list is empty
    else
        xline(x(indices(1)), 'Color', color, 'LineWidth', 1); % Plot the first line for legend
        h = plot(NaN, NaN, 'Color', color, 'LineWidth', 1); % Invisible line for legend
        for index = indices(2:end) % Start from second element to avoid duplicating first line
            xline(x(index), 'Color', color, 'LineWidth', 1);
        end
    end
end
