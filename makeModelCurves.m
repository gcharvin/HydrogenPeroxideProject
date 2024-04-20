

%% steady state concentrations 

% Load or define required functions and parameters
param = get_params(); % Assuming get_params is a function that sets up your model parameters

% Define the keys used for simulation outputs
keys = {'H2O2', 'AntiOxi', 'R_intact', 'P_body', 'R_damaged', 'growth'};

% Initialize a structure to store simulation results for each key
steady_concentrations = struct();
for k = 1:length(keys)
    steady_concentrations.(keys{k}) = [];
end

% Simulation over a range of gamma values
gamma_list = logspace(-1, 1, 32);
for gamma = gamma_list
    param.gamma = gamma;
    [t, steady, converged] = FindAttractor(param); % Assuming this function is correctly translated
    if ~converged
        error('The system did not converge, gamma = %f', gamma);
    end
    for i = 1:length(keys)-1
        steady_concentrations.(keys{i}) = [steady_concentrations.(keys{i}), steady(end, i)];
    end
    steady_concentrations.growth = [steady_concentrations.growth, growth_rate(param.v, steady(end, 3), param.fr_max, param.mu0, param.mR, param.mA)];
end

% Plotting results
figure('Position', [100, 100, 1000, 300]);
for i = 1:length(keys)
    subplot(1, 7, i);
    semilogx(gamma_list, steady_concentrations.(keys{i}));
    set(gca, 'XScale', 'log');
    xlabel('gamma');
    if i == 1
        ylabel('concentration');
    end
    title(keys{i});
end

%% response to step

% Set parameters
param = get_params();  % Assuming a MATLAB function that sets up parameters
param.tramp = 0;
tmax = 200;  % Equivalent to 6 hours

% Prepare to store solutions
sols = {};  % Use cell array to store different solutions

% Simulation values for gamma
gamma_values = [1, 10];

for gamma = gamma_values
    param.gamma = gamma;
    param.Hext_max = 0.0;  % Compute the steady state without H2O2 external
    
    [t, steady, converged] = FindAttractor(param);  % MATLAB function to find attractor

    if ~converged
        error('System did not converge for gamma = %f', gamma);
    end

    % Increase H2O2 external and compute dynamics
    param.Hext_max = 2.0;  % Increase external H2O2
    tlist = logspace(-2, log10(tmax), 512);  % Time points for simulation
    sol = ComputeODE(param, steady, tlist);  % Assuming ComputeODE is a MATLAB function
    sols{end+1} = sol;  % Store solution
end

% Plotting
figure('Position', [100, 100, 1800, 300]);  % Set figure size
keys = {'H2O2', 'AntiOxi', 'R_intact', 'P-body', 'R_damaged'};
for i = 1:length(keys)
    subplot(1, 5, i);
    hold on;
    for j = 1:length(sols)
        plot(sols{j}.x, sols{j}.y(i, :));
    end
    hold off;
    xlabel('Time (min)');
    if i == 1
        ylabel('Concentration');
        legend(arrayfun(@(g) sprintf('gamma = %d', g), gamma_values, 'UniformOutput', false));
    end
    title(keys{i});
end


%%  plot Tdeath 

% Initialization
param = get_params(); % Assume this function sets up your parameters
param.tramp = 0;
gamma_values = [1, 10]; % Gamma values to iterate over
Hext_list = logspace(-2, 3, 256); % External H2O2 concentrations
tlist = logspace(-5, log10(1e5), 4096); % Time points for simulation
Td = cell(1, 2); % Cell array to store T_death for each gamma value

% Simulations
for i = 1:length(gamma_values)
    param.gamma = gamma_values(i);
    param.Hext_max = 0;
    [t, steady, converged] = FindAttractor(param);
    if ~converged
        error('System did not converge for gamma = %f', param.gamma);
    end
    Td{i} = zeros(size(Hext_list)); % Preallocate for speed
    for j = 1:length(Hext_list)
        param.Hext_max = Hext_list(j);
        sol= ComputeODE(param, steady, tlist); % Update to match MATLAB output

        Td{i}(j) = FirstPassageTime(sol, 5, 0.4); % Assume this function is defined
    end
end

% Plotting T_death vs. Hext
figure;
hold on;
for i = 1:length(gamma_values)
    plot(Hext_list, Td{i}, 'DisplayName', sprintf('gamma = %d', gamma_values(i)));
end
plot(logspace(1, 3, 32), 1./logspace(1, 3, 32), 'k--', 'DisplayName', 'slope -1');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('H2O2_{ext}');
ylabel('T_{death}');
legend show;

% Plotting death rate
figure;
hold on;
for i = 1:length(gamma_values)
    validIndices = ~isnan(Td{i}) & Td{i} ~= 0; % Exclude NaNs and zeros
    plot(Hext_list(validIndices), 1./Td{i}(validIndices), 'DisplayName', sprintf('gamma = %d', gamma_values(i)));
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('H2O2_{ext}');
ylabel('Death rate');
legend show;





