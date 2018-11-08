% Function for plotting the three ellipses, given a set of metrics
function plotEllipses(ellipse_properties)

% Define colours to plot the three ellipses
colours = [    [1, 0, 0];   % Red
    [1, 0, 1];   % Purple (magenta)
    [0, 0, 1]    % Blue
    ];

% Convert the vector of all ellipse properties into separate
% ellipse parameters
directions = ellipse_properties([1,4,7]) * pi/180;    % Convert to radians because RegionProps uses degrees
major_lengths = ellipse_properties([2,5,8]);
minor_lengths = ellipse_properties([3,6,9]);

% Plot each ellipse using polar co-ordinates, so vary a parameter,
% t, over [0, 2pi]
t = linspace(0, 2*pi, 500);

% Loop over the three ellipses to plot
for k = 1:3
    
    % Interpret the ellipse in polar co-ordinates, then convert
    % back into (x,y) co-ordinates
    r = 1 ./ sqrt(  ( cos(t - directions(k)) * 2 / major_lengths(k) ).^2  +    ( sin(t-directions(k)) * 2 / minor_lengths(k) ).^2  );
    x = r .* cos(t);
    y = r .* sin(t);
    
    % Plot using (x,y) co-ordinates
    hold on;
    plot(x, y, 'LineWidth', 1.5, 'Color', colours(k,:));
    
    
end
