function visualiseSnapshot(filename, meshfile, vis_mesh)

% Load in the snapshot data
[~, ~, APD, AT, ~] = load_results([pwd,'\Simulations\'], filename);

% Also load in the locations of the points
fileobj_node = fopen([pwd,'/Meshfiles/',meshfile,'.node']);
nodes = textscan(fileobj_node, '%f %f %f %f %f %f','CollectOutput',1);
nodes = nodes{:};
nodes = nodes(2:end, [2 3]);

% Read out the dimensions of the plot
xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));

% Also read out the number of data points in x and y directions
second_row_node = find(nodes(:,2) > ymin, 1);                     % Find the first point above the first row of points
Nx = second_row_node - 2;                                         % The bottom row consists of Nx+1 nodes, so subtract 2 from first point of second row to get the value of Nx
Ny = round( (ymax - ymin) / (nodes(second_row_node,2) - ymin) );  % Ny is found by seeing how many node separations make up the whole y range


% Initialise figure - 
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Create a base set of axes
main_ax = axes('Position', [0.1 0.1 0.8 0.8*Ny/Nx]);

% If visualising the mesh, do it once on a set of overlaid axes
if vis_mesh
    
    % Create axes
    edge_ax = axes('Position', get(main_ax, 'Position'));
    
    % Read out boundary information from file
    fileobj_edge = fopen([pwd,'/Meshfiles/',meshfile,'.edge']);
    edges = textscan(fileobj_edge, '%f %f %f %f','CollectOutput',1);
    edges = edges{:};
    edges = edges(2:end, [2 3]);

    % Convert edges into lines to be plotted
    plot_edges_X = [nodes(edges(:,1), 1), nodes(edges(:,2), 1)];
    plot_edges_Y = [nodes(edges(:,1), 2), nodes(edges(:,2), 2)];
    
    % Plot the edges as red lines
    plot(edge_ax, plot_edges_X', plot_edges_Y', 'r', 'LineWidth', 1.5);
    
    % Remove axis labels and make it transparent, after plotting the data
    set(edge_ax, 'color', 'none', 'xtick', [], 'ytick', [], 'xticklabel', [], 'yticklabel', []);
    % Turn off the axis itself
    %axis(edge_ax, 'off', 'equal');
    % Set the axis limits after plotting the data
    set(edge_ax, 'XLim',[xmin xmax], 'YLim',[ymin ymax]);
    
end

% Visualise APDs successively
for k = 1:size(APD,2)
        
    scatter(main_ax, nodes(:,1), nodes(:,2), 5, APD(:,k));
    colorbar('peer',main_ax, 'EastOutside');
    caxis(main_ax,[50 150]);
    % Set the axis limits after plotting the data
   % axis(main_ax, 'equal');
    set(main_ax, 'XLim',[xmin xmax], 'YLim',[ymin ymax]);
    set(main_ax, 'Position', [0.1 0.1 0.8 0.8*Ny/Nx]);  
    pause(0.25);
    %cla(main_ax);
    
end

% Visualise ATs successively
for k = 1:size(AT,2)
   
    hold on;
    
    scatter(main_ax, nodes(:,1), nodes(:,2), 5, AT(:,k) - min(AT(:,k)));
    caxis(main_ax, [0 30]);
    
    pause(0.25);
    cla(main_ax);
    
end

    
end
