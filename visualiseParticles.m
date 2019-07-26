function visualiseParticles(thetas, theta_mins, theta_maxs, parameter_names, Ds)
% This function visualises a set of particles from a run of an SMC
% algorithm. The code is set up to be used in general, with the user
% supplying necessary information (e.g. parameter names) in the wrapper
% function


% Read out number of parameters, and calculate the number of combinations
N_thetas = size(thetas, 2);
N_combos = N_thetas * (N_thetas - 1) / 2;






%%% PLOT MARGINALS FOR EACH PARAMETER

figure('units','Normalized','OuterPosition',[0 0 1 1]);   % Create a new figure (old ones will be closed by calling function)

% Calculate subplot dimensions based on number of parameters
N_rows = ceil(sqrt(N_thetas / 2));
N_cols = N_rows * 2;

% Loop over each parameter, plotting the marginal distribution for each
for i = 1:N_thetas
   
    % Set up subplot
    subplot(N_rows, N_cols, i);
    
    % Calculate marginal distribution
    [p_theta, theta] = ksdensity(thetas(:,i));
    
    % Plot marginal distribution
    plot(theta, p_theta, 'LineWidth', 2);
    
    % Label plot
    if nargin > 3
        xlabel(parameter_names{i}, 'FontSize', 16);
    else
        xlabel(['\theta_',num2str(i)], 'FontSize', 16);
    end
    ylabel('Probability Density', 'FontSize', 16);
    
    % Define axes
    xlim([theta_mins(i) theta_maxs(i)]);
    
end



%%% PLOT BIVARIATE SCATTERS FOR DIFFERENT PARAMETER COMBINATIONS

figure('units','Normalized','OuterPosition',[0 0 1 1]);    % Create a new figure (old ones will be closed by calling function)

% Calculate subplot dimensions
N_rows = ceil(sqrt(N_combos / 2));
N_cols = N_rows * 2;

% Loop over all combinations
count = 0;
for i = 1:N_thetas-1
    for j = i+1:N_thetas
        
        % Increment counter and select subplot position
        count = count + 1;
        
        % Set up subplot
        subplot( N_rows, N_cols, count);
        
        % Plot particle locations
        scatter(thetas(:,i), thetas(:,j), 10, Ds, 'filled');

        
        % Define axes
        xlim([theta_mins(i), theta_maxs(i)]);
        ylim([theta_mins(j), theta_maxs(j)]);
        
        % Label plot
        if nargin > 3
            xlabel(parameter_names{i}, 'FontSize', 16);
            ylabel(parameter_names{j}, 'FontSize', 16);
        else
            xlabel(['\theta_',num2str(i)], 'FontSize', 16);
            xlabel(['\theta_',num2str(j)], 'FontSize', 16);
        end
        
    end
end


end

