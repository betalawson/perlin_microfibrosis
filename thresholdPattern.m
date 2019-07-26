function presence = thresholdPattern(field, density)

% This function takes an input map of values and converts them to a
% presence/absence configuration according to the requested density. The
% noise pattern is assumed to vary over [0,1], so these are the starting
% limits for the bisection method. A hard limit on the number of steps used
% is implemented to ensure the method does not stall in cases where it is
% impossible to find a threshold that satisfies the specified tolerance
% (due to discrete nature of the problem)

% Define parameters for bisection method
max_iters = 40;
tol = 1e-3;

% Initialise bisection method
window_threshold = [max(field(:)); min(field(:))];

% Loop until iteration limit hit, or solution found within desired
% tolerance
iters = 0;
err = 1e5;
while iters < max_iters && err > tol
    
    % Check midpoint
    check_threshold = ( window_threshold(1) + window_threshold(2) ) / 2;
    check_density = sum(sum(field >= check_threshold)) / numel(field);
    
    % Update window accordingly (recall that as threshold increases,
    % density decreases)
    if check_density < density
        window_threshold(1) = check_threshold;
    else
        window_threshold(2) = check_threshold;
    end
    
    % Calculate error and increment iteration count
    err = abs(density - check_density);
    iters = iters + 1;
    
end

% After loop completes, use the found threshold to convert the noise field
% into a presence/absence pattern
presence = double(field > check_threshold);