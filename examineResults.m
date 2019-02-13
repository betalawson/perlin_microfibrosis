function examineResults(filename)
% This function allows for results to be visualised directly from a .mat
% file. The user specifies the filename of the .mat file to load, and the
% relevant information associated with those results is read directly from
% the filename

% Define parameter ranges - assumes by default
if regexp(filename,'\w*nofibres\w*')
    params_mins = [ 0, 0.01, 0, 1, log(1/2), -pi/2 ];
    params_maxs = [ 0.5, 2, 0.99, 8, log(50), pi/2 ];
else
    params_mins = [ 0, 0.3, 0, 0.01, 0, 1, log(1/2), -pi/2 ];
    params_maxs = [ 0.4, 2, 0.5, 2, 0.99, 8, log(50), pi/2 ];
end

% Load in .mat file data
load([filename,'.mat']);

% Load in images
patterns = FibrosisImageLoader(1);

% Read out the pattern corresponding to this file
if regexp(filename,'\w*interstitial\w*')
    pattern = patterns{1};
elseif regexp(filename,'\w*compact\w*')
    pattern = patterns{2};
elseif regexp(filename,'\w*diffuse\w*')
    pattern = patterns{3};
elseif regexp(filename,'\w*patchy\w*')
    pattern = patterns{4};
elseif regexp(filename,'\w*params\w*')
    load('fibro_seedinfo.mat','permute_tables','offset_tables');
    nameloc = regexp(filename,'\w*params\w*');
    target_density = 0.2;
    points = buildMesh(500,500,1/136);
    switch filename(nameloc+6)     % Reads out character after 'params'
        case '1'
            target_params = [0.15, 0.75, 0.2, 0.7, 0.65, 3, 2.5, -pi/3];
        case '2'
            target_params = [0, 0.75, 0.4, 1.1, 0.9, 3, 3, pi/6];
        case '3'
            target_params = [0.36, 1.1, 0.05, 0.1, 0.3, 5, 1.2, pi/4];
        case '4'
            target_params = [0.3, 0.5, 0.2, 1.8, 0.2, 4, 2, -pi/9];
    end
    pattern = createFibroPattern(points, target_density, target_params, permute_tables{1}, offset_tables{1});  
else
    error('Failed to read out pattern. Check filename for a pattern name.');
end

% Calculate the appropriate ellipse metrics for this pattern
if regexp(filename,'\w*full\w*')
    pattern_M = ellipseMetrics(pattern, [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
elseif regexp(filename,'\w*standard\w*')
    pattern_M = ellipseMetrics(pattern, [0.1 0.5 0.8]);
elseif regexp(filename,'\w*expanded\w*')
    pattern_M = ellipseMetrics(pattern, [0.2 0.6 0.9]);
else
    error('Failed to detect the type of ellipse metrics used. Filename should contain ''FE'', ''SE'' or ''EE''.');
end


% Plot results using the particle data, and the target data
visualiseResults_Fibrosis(particles.thetas, particles.vals, particles.metrics, particles.Ds, pattern, pattern_M,params_mins,params_maxs);