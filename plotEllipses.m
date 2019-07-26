% Function for plotting the three ellipses, given a set of metrics
function plotEllipses(ellipse_properties)

% Define colours to plot the three ellipses
% colours = [    
%     [0.85, 0, 0];            % Dark red
%     [1, 0.5, 0.5];           % Pink
%     [0.5, 0.5, 1];           % Light blue
%     [0, 0, 0.85];            % Dark blue
%     [0, 0, 0];               % Black
%     [0.3, 0.3, 0.3];         % Dark grey
%     [0.75, 0.75, 0.75];      % Light grey
%     [0.6, 0.6, 0.3];         % Dark olive
%     [0.9, 0.9, 0.6];         % Light olive
%     ];
colours = [    
    [0.85, 0, 0];            % Dark red
    [0.9, 0.25, 0.25];       % Moderate-dark red
    [0.95, 0.5, 0.5];        % Moderate red
    [1, 0.75, 0.75];         % Light red
    [0.85, 0.75, 0.85];      % Light purple
    [0.75, 0.75, 1];         % Light blue
    [0.5, 0.5, 0.95];        % Moderate blue
    [0.25, 0.25, 0.9];       % Moderate-dark blue
    [0, 0, 0.85];            % Dark blue
    ];

% Convert the vector of all ellipse properties into separate
% ellipse parameters
directions = ellipse_properties(1:3:length(ellipse_properties)) * pi/180;    % Convert to radians because RegionProps uses degrees
major_lengths = ellipse_properties(2:3:length(ellipse_properties));
minor_lengths = ellipse_properties(3:3:length(ellipse_properties));

% Plot each ellipse using polar co-ordinates, so vary a parameter,
% t, over [0, 2pi]
t = linspace(0, 2*pi, 500);

% Loop over the three ellipses to plot
for k = 1:length(ellipse_properties)/3
    
    % Interpret the ellipse in polar co-ordinates, then convert
    % back into (x,y) co-ordinates
    r = 1 ./ sqrt(  ( cos(t - directions(k)) * 2 / major_lengths(k) ).^2  +    ( sin(t-directions(k)) * 2 / minor_lengths(k) ).^2  );
    x = r .* cos(t);
    y = r .* sin(t);
    
    % Plot using (x,y) co-ordinates
    hold on;
    plot(x, y, 'LineWidth', 1.5, 'Color', colours(k,:));
    
end
