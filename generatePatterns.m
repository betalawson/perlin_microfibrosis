function patterns = generatePatterns(params, density, N_patterns, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
% (mesh) - optionally provided mesh to specify size of patterns
%
% PARAMETER INFORMATION:
%
% 1 - FIBRENESS: The extent to which patterns exhibit long, thin fibres 
%     aligned in consistent directions
%     ::: If set to NaN, a pattern without fibres will be created :::
% 2 - FIBRE SEPARATION: The average spacing between fibres (in units
%     matching those used in input mesh
% 3 - PATCHINESS: The extent of inconsistency in pattern density (high
%     patchiness will produce distinct regions of higher and lower density)
% 4 - FEATURE SIZE: The overall size of obstacle features in the mesh (in
%     units matching the mesh)
% 5 - ROUGHNESS: The roughness of feature edges (values from [0,1], may
%     cause issues for values of precisely 0 or 1)
% 6 - PATCH SIZE: The size of regions over which density varies (with 
%     extent according to PATCHINESS)
% 7 - FIBRE ALIGNMENT: The extent to which non-fibre features are aligned
%     to the fibre direction (i.e. extent of feature anisotropy)
% 8 - DIRECTION: An angle (in radians) defining the orientation of fibres
%     and/or feature anisotropy

% Load in the seed data
load('fibro_seedinfo.mat', 'permute_tables', 'offset_tables');

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 4
    mesh = buildMesh(250, 400, 1/136);
end

% Create the requested number of patterns
for m = 1:N_patterns
    
    % Use the fibre-free generator if NaNs are present in input params
    % vector, or if only non-fibre parameters provided, otherwise 
    % use the standard generator
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, params(3:8), permute_tables{m}, offset_tables{m});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, params, permute_tables{m}, offset_tables{m});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, density, params, permute_tables{m}, offset_tables{m});
    end
    
    % Store this pattern
    patterns{m} = presence;
    
end


% Initialise a figure to plot some of the patterns
figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Plot the first ten patterns (or less if less were requested)
for k = 1:min(10, N_patterns)
   
    % Plot this pattern
    subplot(2,5,k);
    imagesc(patterns{k}); axis('equal', 'off');
    colormap(fibroclr);
    
end