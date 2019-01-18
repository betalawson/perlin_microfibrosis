function [presence, S, P_f, P_p] = createFibroPattern(mesh, density, params, Ps, offsets)

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
[fibreness, fibre_sep, patchiness, feature_size, roughness, patch_size, fibre_alignment, direction] = deal(params{:});


% Create a rotated set of points for the application of anisotropy and
% creation of fibre-aligned pattern. Stored as two rows for ease of matrix
% transforms and input into C functions
R_points = [ [ cos(direction) sin(direction) ]; [-sin(direction), cos(direction)] ] * mesh.points';


% Create new permutation tables from the provided by applying it to itself,
% and then again
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k,Ps(k,:)+1);
    Ps3(k,:) = Ps(k,Ps2(k,:)+1);
end


%%% SINUSOIDAL NOISE FOR FIBROUS PATTERNS

% Define parameters for the phase field (decided via experimentation)
n_fibres_similarity = 0.8;
wiggle_feature_length = 0.8;
phasefield_strength = 1;

% Create an octave noise field with these properties (first transform
% points, then call Octave2D)
%phasefield_points = [ R_points(1,:) / sqrt(phasefield_anisotropy); R_points(2,:) * sqrt(phasefield_anisotropy) ];      % Scale dimensions according to anisotropy (multiplication by "D" in paper)
%phasefield_points = phasefield_points / phasefield_featuresize;                                                                % Attain base featuersize for this noisefield (appliation of 1/d factor in paper)
phasefield_points = [ R_points(1,:) / wiggle_feature_length; R_points(2,:) / (n_fibres_similarity * fibre_sep) ];      % Scale dimensions according to anisotropy (multiplication by "D" in paper)
phasefield = Octave2D(phasefield_points, 4, 0.5, Ps, offsets);    % (4 octaves, low roughness for phasefield)
% Create the sinusoidal pattern using cos(RX^T [0; 1] ), then modulate phase using
% the created phasefield
S = 0.5 + 0.5 * cos(2*pi * (R_points(2,:) / fibre_sep + phasefield_strength * (phasefield - 0.5) ) );

% Sharpen this field by powering it up many times
S = S .* S .* S .* S .* S .* S .* S .* S .* S .* S .* S .* S .* S .* S .* S;



%%% CREATE THE MAIN FIBROSIS DEPOSIT EFFECT

% Transform points according to input parameters, then call Octave2D
P_f_points = [ R_points(1,:) / sqrt(fibre_alignment); R_points(2,:) * sqrt(fibre_alignment) ];
P_f = Octave2D( P_f_points / feature_size, 4, roughness, Ps2, offsets);



%%% CREATE A LARGE-SCALE PERLIN NOISE PATTERN FOR DENSITY VARIATION

% Use Octave2D with scaling of point co-ords to attain desired patch_size
P_p = Octave2D(mesh.points' / patch_size, 4, 0.5, Ps3, offsets);



%%% TAKE A COMBINATION OF THESE NOISEFIELDS TO GET THE FINAL PATTERN
%noise = (1 - patchiness) * ( fibreness * S + (1 - fibreness) * P_f ) + patchiness * P_p;
noise = ( 1 - fibreness + fibreness * S ) .* P_f + patchiness * P_p;

%%% THRESHOLD THIS NOISE TO GET PRESENCE/ABSENCE OF REQUESTED DENSITY
presence = thresholdPattern(noise, density);


%%% CONVERT BACK TO MATRICES
presence = reshape(presence', mesh.Nx, mesh.Ny)';
S = reshape(S', mesh.Nx, mesh.Ny)';
P_f = reshape(P_f', mesh.Nx, mesh.Ny)';
P_p = reshape(P_p', mesh.Nx, mesh.Ny)';

end

