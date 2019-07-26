function visualiseResults_Fibrosis(params, patterns, ellipses, discrepancies, target_pattern, target_ellipses, params_mins, params_maxs)
% This function plots the particles in terms of their positions in the
% theta space, the ellipses used for pattern representation, and the
% patterns themselves.

%%% PARAMETER DEFINITION

N_vis = 11;
plot_type = 'random';           % 'best' - picks best particles in terms of discrepancy
                                % 'random' - picks a random set of particles each time the function is called



%%% INITIAL SETUP

% Close all current figures, so that this code can create new ones (allows
% for live visualisation as the SMC-ABC algorithm works)
close all;

% Read out number of particles provided and ensure it's enough to plot the
% requested amount
N = length(patterns);
if N < N_vis
    error('Not enough particles to plot the requested number. Adjust N_vis in plotting code');
end

% Select which particles to use for plotting according to the plotting type
% specified above
switch plot_type
    
    case 'random'    % Select N_vis patterns at random
        
        N = length(patterns);
        R = randperm(N);
        plot_patterns = patterns(R(1:N_vis));
        plot_ellipses = ellipses(R(1:N_vis),:);
        plot_discrepancies = discrepancies(R(1:N_vis));
        
    case 'best'      % Select the best N_vis patterns according to discrepancy
        
        [plot_discrepancies, I] = sort(discrepancies);
        plot_patterns = patterns(I(1:N_vis));
        plot_ellipses = ellipses(I(1:N_vis),:);
        
end

% Calculate the number of rows and columns to use for subplots in
% visualising the patterns and the ellipses
N_rows = ceil(sqrt(N_vis / 2));
N_cols = N_rows * 2;



%%% PLOT PARTICLES USING GENERALISED CODE

% Define parameter names for the fibrosis problem here
if size(params,2) == 8
    parameter_names = {'Fibreness', 'Fibre Separation', 'Patchiness', 'Feature Size', 'Roughness', 'Patch Size', '(log) Anisotropy Ratio', 'Direction'};
elseif size(params,2) == 6
    parameter_names = {'Patchiness', 'Feature Size', 'Roughness', 'Patch Size', '(log) Anisotropy Ratio', 'Direction'};
else
    fprintf('WARNING: Parameter names are likely incorrect!\n');
    parameter_names = {'Fibreness', 'Fibre Separation', 'Patchiness', 'Feature Size', 'Roughness', 'Patch Size', '(log) Anisotropy Ratio', 'Direction'};
end
% Call particle plotter
visualiseParticles(params, params_mins, params_maxs, parameter_names, discrepancies);


%%% CALL PLOTTING SUBFUNCTIONS FOR PATTERNS AND ELLIPSES
figure('units','Normalized','OuterPosition',[0 0 1 1]);
visualisePatterns;
figure('units','Normalized','OuterPosition',[0 0 1 1]);
visualiseEllipses;


    % Function for visualising the patterns themselves. The target pattern
    % is plotted first, then a selection of the patterns associated with
    % current particles
    function visualisePatterns
        
        % Plot the target pattern in the top left, and mark it as such
        subplot(N_rows, N_cols, 1);
        imagesc(target_pattern);
        axis equal;
        title('TARGET', 'FontSize', 20);
        
        % Plot the remaining patterns
        for k = 1:N_vis
            
            subplot(N_rows, N_cols, k+1);
            imagesc(plot_patterns{k});
            axis equal;
            xlabel(['Discrepancy: ',num2str(plot_discrepancies(k))], 'FontSize', 16);
            
        end
        
    end


    % Function for visualising the ellipses associated with each pattern.
    % Ellipses for the target pattern are plotted first, then a selection
    % of the ellipses of the particles
    function visualiseEllipses
        
        % Determine appropriate plotting limits, so each subplot can use
        % consistent axes. The major axis length of the largest (containing
        % most power of spectrum) ellipse is used to define the length of
        % these axes.
        max_major = max( [target_ellipses(end-1); plot_ellipses(1:N_vis,end-1)] );
        ax_start = -max_major * 0.6;
        ax_end = max_major * 0.6;
        
        % Plot the target pattern in the top left, and mark it as such
        subplot(N_rows, N_cols, 1);
        hold on;
        plotEllipses(target_ellipses);
        xlim([ax_start ax_end]);
        ylim([ax_start ax_end]);
        axis equal;
        title('TARGET', 'FontSize', 20);
        
        % Plot the remaining patterns
        for k = 1:N_vis
            
            subplot(N_rows, N_cols, k+1);
            hold on;
            plotEllipses(plot_ellipses(k,:));
            xlim([ax_start ax_end]);
            ylim([ax_start ax_end]);
            axis equal;
            xlabel(['Discrepancy: ',num2str(plot_discrepancies(k))], 'FontSize', 16);
            
        end
        
    end

end

