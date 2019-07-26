function patterns = generatePatterns(params, density, N_patterns, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
% (mesh) - optionally provided mesh to specify size of patterns

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
    % vector, otherwise standard generator
    if any(isnan(params))
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