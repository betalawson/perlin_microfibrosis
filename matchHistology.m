function matchHistology(N_parts, histo_pattern, ellipse_mode, visualise, options)
% This function creates a population of Perlin noise-derived patterns that
% all fall within a sufficient degree of similarity to a target
% histological pattern.
%
% INPUTS:
%
% N_parts:        Number of particles to use in the SMC-ABC method
%
% histo_pattern:  Name of histological pattern to match, options are
%                      'interstitial', 'compact', 'diffuse', 'patchy'
% 
% ellipse_mode:   The set of ellipses to use for metric calculation,
%                 options are
%                       'standard', 'expanded', 'full'                 
%
% visualise:      Flag (0 or 1) for visualisation of SMC-ABC algorithm as
%                 it progresses (visualisation function is specified)
%
% (options):      Options for the SMC-ABC algorithm (see SMC-ABC code for
%                 definitions and default values)

% Define the "pixel width" from the histological data
pixel_width = 1/136;      % (mm)

% Specify the histological pattern to match (if doing so), and whether to use full or reduced patterns
reduce_patterns = 1;                % 0 - take the full images,   1 - reduce to representative sections

% PARAMETERS ARE:
% [ fibreness, fibre_separation, patchiness, feature_size, roughness, patch_size, alignment_ratio]

% Define the prior for ABC-SMC. Mins and maxs are the ranges of prior,
% scale_param selects between uniform and logarithmic prior
params_mins = [ 0, 0.3, 0, 0.01, 0, 1, 1/2, -pi/2 ];
params_maxs = [ 0.4, 2, 0.5, 2, 0.99, 8, 50, pi/2 ];
scale_param = logical([ 0, 0, 0, 0, 0, 0, 1, 0]);


%%% CHECK IF OPTIONS SUPPLIED
if nargin > 4
    options_supplied = 1;
else
    options_supplied = 0;
end

%%% RUN MAIN CODE

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
    otherwise
        error('No histological pattern with that name!');
end
        
% Set up the mesh according to the pattern selected
[Ny, Nx] = size(target_pattern);
points = buildMesh(Nx, Ny, pixel_width);
    
% Calculate density of pattern to be matched
density = sum(target_pattern(:)) / numel(target_pattern);
    

% Calculate limits on the parameters after transforming scale parameters to
% uniform
theta_mins = params_mins;
theta_mins(scale_param) = log(theta_mins(scale_param));
theta_maxs = params_maxs;
theta_maxs(scale_param) = log(theta_maxs(scale_param));

% Define the function used for simulating the model (here generating a pattern) in the SMC-ABC
f_simulate = @(params, P, offsets) createFibroPattern(points, density, params, P, offsets);

% Define the function used for calculating summary statistics in the SMC-ABC
switch ellipse_mode
    case 'standard'
        f_summaries = @(pattern) ellipseMetrics(pattern, [0.1 0.5 0.8]);
    case 'expanded'
        f_summaries = @(pattern) ellipseMetrics(pattern, [0.2 0.6 0.9]);
    case 'full'
        f_summaries = @(pattern) ellipseMetrics(pattern, [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
    otherwise
        error('Ellipse_mode must be ''standard'', ''expanded'', or ''full''!');
end

% Define the function used for calculating discrepancies between summary statistics
f_discrepancy = @ellipseDiscrepancy;

% Define the function for visualisation
target_ellipses = f_summaries(target_pattern);
f_visualise = @(params, patterns, ellipses, discrepancies) visualiseResults_Fibrosis(params, patterns, ellipses, discrepancies, target_pattern, target_ellipses, theta_mins, theta_maxs);

% Run the SMC-ABC
if options_supplied
    [part_thetas, part_vals, part_metrics, part_Ds] = performABCSMC_Fibro(N_parts, f_simulate, f_summaries, f_discrepancy, target_pattern, params_mins, params_maxs, scale_param, 1e-2, visualise, f_visualise, options);
else
    [part_thetas, part_vals, part_metrics, part_Ds] = performABCSMC_Fibro(N_parts, f_simulate, f_summaries, f_discrepancy, target_pattern, params_mins, params_maxs, scale_param, 1e-2, visualise, f_visualise);
end

% Save data to a struct
particles.thetas = part_thetas;
particles.vals = part_vals;
particles.metrics = part_metrics;
particles.Ds = part_Ds;

% Set up a filename according to options supplied
filename = [histo_pattern, num2str(N_parts),'_', ellipse_mode];

% Append options to filename if supplied
if options_supplied
        
    % Check all fields and append to filename
    if isfield(options,'jumping_type')
        filename = [filename,'_jt_',options.jumping_type];
    end
    if isfield(options,'resample_weighting')
        filename = [filename,'_rw_',options.resample_weighting];
    end
    if isfield(options,'metric_weighting')
        filename = [filename,'_mw_',options.metric_weighting];
    end
    if isfield(options,'keep_fraction')
        filename = [filename,'_kf_',num2str(options.keepfraction)];
    end
    if isfield(options,'max_MCMC_steps')
        filename = [filename,'_maxMCMC_',num2str(options.max_MCMC_steps)];
    end
    
end

% Save the particles
save([filename,'.mat'],'particles');


end

