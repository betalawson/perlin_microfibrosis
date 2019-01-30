function visualiseSim(filename, meshfile, vis_mesh, save_anim)
% This function visualises the simulation results contained in the full
% data .h5 files produced by Chaste. The meshfile is also provided, so that
% the underlying mesh can be visualised.

% Load the results with the given filename - looks for simulations in a
% /Simulations/ folder
[t, V, ~, ~, ~] = load_results([pwd,'\Simulations\'], filename);

meshfile = [pwd,'\Meshfiles\',meshfile];

% Read out the node information from the .node meshfile
fileobj_node = fopen([meshfile,'.node']);
nodes = textscan(fileobj_node, '%f %f %f %f %f %f','CollectOutput',1);
nodes = nodes{:};
nodes = nodes(2:end, [2 3]);

xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));

% If plotting mesh information, read out this information from the element
% and edge files
if vis_mesh

    % Read out element information from the .ele meshfile
    fileobj_ele = fopen([meshfile,'.ele']);
    elements = textscan(fileobj_ele, '%f %f %f %f','CollectOutput',1);
    elements = elements{:};
    elements = elements(2:end, [2 3 4]);

    % Convert element information into triangles to be plotted
    plot_tri_X = [nodes(elements(:,1), 1), nodes(elements(:,2), 1), nodes(elements(:,3), 1) ];
    plot_tri_Y = [nodes(elements(:,1), 2), nodes(elements(:,2), 2), nodes(elements(:,3), 2) ];

    % Read out boundary information from file
    fileobj_edge = fopen([meshfile,'.edge']);
    edges = textscan(fileobj_edge, '%f %f %f %f','CollectOutput',1);
    edges = edges{:};
    edges = edges(2:end, [2 3]);

    % Convert edges into lines to be plotted
    plot_edges_X = [nodes(edges(:,1), 1), nodes(edges(:,2), 1)];
    plot_edges_Y = [nodes(edges(:,1), 2), nodes(edges(:,2), 2)];
    
end

% If saving an animation, set it up here
if save_anim
    Vidobj = VideoWriter([filename, '.avi']);
    Vidobj.FrameRate = 15;
    open(Vidobj);
end

% Create a basic colormap that displays unexcited cells as dark blue, and
% excitement as red, through a grey
fibro_colormap = [ [ linspace(0,0.5,81)', linspace(0,0.5,81)', linspace(0.4,0.7,81)' ];
                   [ linspace(0.5,1,40)', linspace(0.5,1,40)', linspace(0.7,0.5,40)' ] ]; 
                  
                   
% Create a new figure to display the data
figure('units','normalized','OuterPosition',[0 0 1 1]);

% Close open files (just avoiding fclose all for better compatibility
% outside this code)
fclose(fileobj_node);
if vis_mesh
    fclose(fileobj_ele);
    fclose(fileobj_edge);
end

% Loop over timesteps
for i = 1:size(V,2)
    
    % Visualise data
    hold on;
    
   
    % Visualise the nodes themselves, colourised according to voltage value
    scatter(nodes(:,1), nodes(:,2), 5, V(:,i), 'filled');
    caxis([-100 20]);
    axis equal;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    %colormap(fibro_colormap);
    
    % Visualise the elements and boundaries if requested
    if vis_mesh
    
        % Visualise element triangles
        %patch(plot_tri_X', plot_tri_Y', 90);
    
        % Visualise boundaries
        plot(plot_edges_X', plot_edges_Y', 'r', 'LineWidth', 1.5);
        
    end
    
    title(['Time t = ',num2str(t(i)),'ms']);
    
    % Update frame, or save it
    if ~save_anim
        pause(0.5);
    else
        frame = getframe(gcf);
        writeVideo(Vidobj,frame);
    end
    
    % Clear the image because we're using hold on, except final frame
    if i ~= size(V,2)
        cla; clf;
    end
    
end

% Close video object if one was created
if save_anim
    close(Vidobj);
end
