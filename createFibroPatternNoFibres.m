function [presence, O_b, O_d] = createFibroPatternNoFibres(mesh, density, params, Ps, offsets)
% This function takes a list of points and a set of parameters (listed 
% below), and creates a pattern of fibrosis accordingly.
%
% Usage:    presence = CreateFibroPattern(points, direction, density, params, Ps, offsets)
%
% INPUTS:   mesh:      information that defines the mesh over which to create a pattern
%           density:   the density of fibrosis
%           params:    a set of parameters defining the noise pattern (see below)
%           Ps:        an m x n matrix, which is m rows of the numbers 0:n-1 arranged in random order (permutation tables for random assignment of vectors in Perlin noise)
%           offsets:   an m x 2 matrix that specifies grid offsets for each octave in octave noise
%
%           ( m is the maximum number of octaves that will be requested by
%           this function, and n is the number of 
%
% OUTPUTS:  presence:  a presence/absence map of fibrosis
%           (S):       the sinusoidal field used for fibre-type patterns (optional)
%           (Pf):      the Perlin noise field for fibrosis (optional)
%           (Pp):      the Perlin noise field for density variation (optional)
%
% PARAMS:   Parameters are provided as a single vector, for ease of
%           interface with other code.
%           Params vector is specified as:
%           [ fibreness, fibre_separation, patchiness, feature_size, roughness, patch_size, alignment_ratio, direction]

% Read out paramaters from vector
params = num2cell(params);
[patchiness, feature_size, roughness, patch_size, fibre_alignment, direction] = deal(params{:});


% Create a rotated set of points for the application of anisotropy and
% creation of fibre-aligned pattern. Stored as two rows for ease of matrix
% transforms and input into C functions
R_points = [ [ cos(direction) sin(direction) ]; [-sin(direction), cos(direction)] ] * mesh.points';


% Create new permutation tables from the provided by applying it to itself
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k,Ps(k,:)+1);
end


%%% CREATE THE MAIN FIBROSIS DEPOSIT EFFECT

% Transform points according to input parameters, then call Octave2D
P_f_points = [ R_points(1,:) / sqrt(fibre_alignment); R_points(2,:) * sqrt(fibre_alignment) ];
O_b = Octave2D( P_f_points / feature_size, 4, roughness, Ps, offsets);



%%% CREATE A LARGE-SCALE PERLIN NOISE PATTERN FOR DENSITY VARIATION

% Use Octave2D with scaling of point co-ords to attain desired patch_size
O_d = Octave2D(mesh.points' / patch_size, 3, 0.5, Ps2, offsets);



%%% TAKE A COMBINATION OF THESE NOISEFIELDS TO GET THE FINAL PATTERN
noise = O_b + patchiness * O_d;

%%% THRESHOLD THIS NOISE TO GET PRESENCE/ABSENCE OF REQUESTED DENSITY
presence = thresholdPattern(noise, density);


%%% CONVERT BACK TO MATRICES
presence = reshape(presence', mesh.Nx, mesh.Ny)';
O_b = reshape(O_b', mesh.Nx, mesh.Ny)';
O_d = reshape(O_d', mesh.Nx, mesh.Ny)';

%%% FLIP TO CONVERT BACK TO IMAGES (starting at top left instead of bottom
%%% left)
presence = flipud(presence);
O_b = flipud(O_b);
O_d = flipud(O_d);

end
