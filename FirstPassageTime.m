function t_passage = FirstPassageTime(sol, variable_index, threshold)
    % Extract the time and variable data from sol
    t = sol.x;
    variable_data = sol.y(variable_index, :);

    % Find the first index where the variable data exceeds the threshold
    idx = find(variable_data > threshold, 1, 'first');

    % Check if such an index exists
    if isempty(idx)
        t_passage = NaN; % Return NaN if the threshold is never crossed
    else
        t_passage = t(idx); % Return the corresponding time value
    end
end
