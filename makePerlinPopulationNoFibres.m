function makePerlinPopulationNoFibres(filename)
% This function creates a population of Perlin noise-derived patterns that
% all fall within a sufficient degree of similarity to a target pattern,
% which may be one of the histological 

% Define the "pixel width", taken from the histological images but also
% used to give physical dimensions to pixels in synthetic data
pixel_width = 1/136;      % (mm)

% Specify the number of particles to use for ABC-SMC
N_parts = 250;

% Specify whether to match a histological pattern, or to match to a set of provided "true" parameter values
match_histo = 1;
% Specify the histological pattern to match (if doing so), and whether to use full or reduced patterns
histo_pattern = 'diffuse';     % Options are 'interstitial', 'compact', 'diffuse', 'patchy'
reduce_patterns = 1;                % 0 - take the full images,   1 - reduce to representative sections

% Specify the properties of the target pattern, used when match_histo is set to false
Nx = 400;
Ny = 400;
density = 0.2;



% PARAMETERS ARE:
% [ patchiness, feature_size, roughness, patch_size, alignment_ratio, direction]

% Define the prior for ABC-SMC. Mins and maxs are the ranges of prior,
% scale_param selects between uniform and logarithmic prior
params_mins = [ 0, 0.1, 0, 1, 1/2, -pi/2 ];
params_maxs = [ 0.5, 2, 1, 8, 50, pi/2 ];
scale_param = logical([ 0, 0, 0, 0, 1, 0]);

% Parameters to match when using ABC-SMC to re-match a provided pattern
match_params = [ 0.2, 0.75, 0.85, 2.8, 3, -5*pi/8];



%%% RUN MAIN CODE

if match_histo   % Attempting to fit to one of the four histological patterns

	% Load in the requested pattern
	presences = FibrosisImageLoader(reduce_patterns);
    switch histo_pattern
        case 'interstitial'
            target_pattern = presences{1};
		case 'compact'
            target_pattern = presences{2};
		case 'diffuse'
            target_pattern = presences{3};
		case 'patchy'
            target_pattern = presences{4};	
    end
        
    % Set up the mesh according to the pattern selected
    [Ny, Nx] = size(target_pattern);
    points = buildMesh(Nx, Ny, pixel_width);
    
    % Calculate density of pattern to be matched
    density = sum(target_pattern(:)) / numel(target_pattern);
    

else   % Testing ability to re-capture given parameters

    % Set up the mesh according to the properties specified
    points = buildMesh(Nx, Ny, pixel_width);
    
    % Create the pattern using the set of parameters to match
    target_pattern = createFibroPattern(points, density, match_params);
    
end

% Calculate limits on the parameters after transforming scale parameters to
% uniform
theta_mins = params_mins;
theta_mins(scale_param) = log(theta_mins(scale_param));
theta_maxs = params_maxs;
theta_maxs(scale_param) = log(theta_maxs(scale_param));

% Define the function used for simulating the model (here generating a pattern) in the ABC-SMC
f_simulate = @(params, P, offsets) createFibroPatternSeededNoFibres(points, density, params, P, offsets);
% Define the function used for calculating summary statistics in the ABC-SMC
f_summaries = @ellipseMetrics;
% Define the function for visualisation
target_ellipses = f_summaries(target_pattern);
f_visualise = @(params, patterns, ellipses, discrepancies) visualiseResults_Fibrosis(params, patterns, ellipses, discrepancies, target_pattern, target_ellipses, theta_mins, theta_maxs);

% Run the ABC-SMC
%[part_thetas, part_vals, part_metrics, part_Ds] = performABCSMC(N_parts, f_simulate, f_summaries, target_pattern, params_mins, params_maxs, scale_param, 1e-2, 1, f_visualise);
[part_thetas, part_vals, part_metrics, part_Ds] = performABCSMC_FibroSeeds(N_parts, f_simulate, f_summaries, target_pattern, params_mins, params_maxs, scale_param, 1e-2, 1, f_visualise);

% (TO CHANGE TO MORE BAYESIAN INTERPRETATION)
% Save only the unique particles
[particles.thetas, I] = unique(part_thetas, 'rows');
particles.vals = part_vals(I);
particles.metrics = part_metrics(I,:);
particles.Ds = part_Ds(I);

% Check if a filename was provided, otherwise default it to 'particles'
if nargin == 0
    filename = 'particles';
end

% Save the particles
save([filename,'.mat'],'particles');


end

