function [presence, P_f, P_p] = createFibroPatternSeededNoFibres(mesh, density, params, P, offsets)

% This function takes a list of points and a set of parameters (listed 
% below), and creates a pattern of fibrosis accordingly.
%
% Usage:    presence = CreateFibroPatternNew(points, direction, density, params, seed)
%
% INPUTS:   mesh:      information that defines the mesh over which to create a pattern
%           density:   the density of fibrosis
%           params:    a set of parameters defining the noise pattern (see below)
%
% OUTPUTS:  presence:  a presence/absence map of fibrosis
%           (Pf):      the Perlin noise field for fibrosis (optional)
%           (Pp):      the Perlin noise field for density variation (optional)
%
% PARAMS:   Parameters are provided as a single vector, for ease of
%           interface with other code.
%           Params vector is specified as:
%           [ patchiness, feature_size, roughness, patch_size, alignment_ratio, direction]

% Read out paramaters from vector
params = num2cell(params);
[patchiness, feature_size, roughness, patch_size, fibre_alignment, direction] = deal(params{:});


% Create a rotated set of points for the application of anisotropy and
% creation of fibre-aligned pattern. Stored as two rows for ease of matrix
% transforms and input into C functions
R_points = [ [ cos(direction) sin(direction) ]; [-sin(direction), cos(direction)] ] * mesh.points';


% Create new permutation tables from the provided by applying it to itself,
% and then again
P2 = P(P+1);


%%% CREATE THE MAIN FIBROSIS DEPOSIT EFFECT

% Transform points according to input parameters, then call Octave2D
P_f_points = [ R_points(1,:) / sqrt(fibre_alignment); R_points(2,:) * sqrt(fibre_alignment) ];
P_f = Octave2DSeeded( P_f_points / feature_size, 4, roughness, P, offsets);



%%% CREATE A LARGE-SCALE PERLIN NOISE PATTERN FOR DENSITY VARIATION

% Use Octave2D with scaling of point co-ords to attain desired patch_size
P_p = Octave2DSeeded(mesh.points' / patch_size, 3, 0.5, P2, offsets);



%%% TAKE A COMBINATION OF THESE NOISEFIELDS TO GET THE FINAL PATTERN
%noise = (1 - patchiness) * ( fibreness * S + (1 - fibreness) * P_f ) + patchiness * P_p;
noise = P_f + patchiness * P_p;

%%% THRESHOLD THIS NOISE TO GET PRESENCE/ABSENCE OF REQUESTED DENSITY
presence = thresholdPattern(noise, density);


%%% CONVERT BACK TO MATRICES
presence = reshape(presence', mesh.Nx, mesh.Ny)';
P_f = reshape(P_f', mesh.Nx, mesh.Ny)';
P_p = reshape(P_p', mesh.Nx, mesh.Ny)';

end

