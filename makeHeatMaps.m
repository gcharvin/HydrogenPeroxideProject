BIN = 16; % Number of bins
param = get_params(); % Assuming a function to set initial parameters
param.tramp = 0;
tmax = 360; % 6 hours

gamma_range = logspace(-0.3, 0.7, BIN);
Hext_range = logspace(-0.7, 0.7, BIN);

% Preallocate the output array for efficiency
% Dimensions: [length of gamma_range, length of Hext_range, 7 outputs per combination]
output = zeros(length(gamma_range), length(Hext_range), 7);

for i = 1:length(gamma_range)
    gamma = gamma_range(i);
    param.gamma = gamma;
    param.Hext_max = 0.0;
    [t, steady, converged] = FindAttractor(param);

    for j = 1:length(Hext_range)
        Hext = Hext_range(j);
        param.Hext_max = Hext;
        tlist = logspace(-2, log10(tmax), 512);
        sol = ComputeODE(param, steady, tlist);
        [~, minIndex] = min(abs(sol.y(5, :) - 0.4));
        Tdeath = sol.x(minIndex);
        
        % Store results directly in the output 3D array
        output(i, j, 1) = gamma;
        output(i, j, 2) = Hext;

        output(i, j, 3:7) = sol.y(:, end).'; % Transpose to make row vector
        output(i, j, 8) = Tdeath;


    end
end

data=output;

% Extract columns
gamma = data(:,:,1); % Gamma values
Hext = data(:,:,2); % External H2O2 concentrations
H_internal = data(:,:,3); % Internal H2O2
Antioxidant = data(:,:,4); % Antioxidant
Ribosome = data(:,:,5); % Intact free Ribosome
RibosomeTotal = data(:,:,5) + data(:,:,6); % Total intact Ribosome
Pbody = data(:,:,6); % P-body
Rdamaged = data(:,:,7); % Damaged Ribosomes
Tdeath = data(:,:,8); % Tdeath, handling cases where Tdeath is not defined (should be handled if data file has any NaN or Inf)

% Set up the figure and axes, use unique figure for each to save separately
vars = {H_internal, Antioxidant, Ribosome, RibosomeTotal, Pbody, Rdamaged, Tdeath};
titles = {'H internal', 'Antioxidant', 'Intact free Ribosome', ...
    'Total intact Ribosome', 'P-body', 'R damaged', 'Tdeath'};
filenames = {'H', 'A', 'R', 'Rtot', 'Pbody', 'Rdamaged', 'Tdeath'};

figure; % Create a single figure
numVars = length(vars); % Number of variables to plot
numCols = ceil(sqrt(numVars)); % Number of columns in subplot grid (adjust layout as needed)
numRows = ceil(numVars / numCols); % Number of rows in subplot grid

for k = 1:numVars
    subplot(numRows, numCols, k); % Position each heatmap in the subplot grid
    set(gca, 'XScale', 'log', 'YScale', 'log'); % Log-log scale
    colormap('hot'); % Use the same colormap for all plots for consistency

    % Handling data for pcolor
    unique_gamma = unique(gamma);
    unique_Hext = unique(Hext);
    [X, Y] = meshgrid(unique_gamma, unique_Hext);

    % Check the number of elements
    expectedNumElements = length(unique_Hext) * length(unique_gamma);

    if numel(vars{k}) == expectedNumElements
        Z = reshape(vars{k}, length(unique_Hext), length(unique_gamma));
        pcolor(X, Y, Z);
        shading interp; % This makes the heatmap look continuous
        colorbar;
        xlabel('{\gamma} (1/min)');
        ylabel('H_{ext} (mM)');
        title(titles{k});
        set(gca, 'ColorScale', 'log'); % Set color axis to logarithmic if needed
    else
        error('The number of elements in vars{%d} does not match expected dimensions.', k);
    end
end

% Optionally, adjust subplot margins
set(gcf, 'Position', [100, 100, 1400, 900]); % Resize the figure to fit all subplots

% Save the entire figure as EPS
%print('figures/AllVars.eps', '-depsc');
