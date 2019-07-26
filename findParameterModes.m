function param_modes = findParameterModes(particles, present_params, scale_params)
% This function takes a particle set, and optional additional arguments
% specifying some information about the data, and outputs 
%
% INPUTS:
%
% particles - a struct type variable containing the set of particles
% (present_params) - optional argument specifying which variables are present in the particle set
% (scale_params) - optional argument specifying which variables are scaling parameters

% Read out the number of parameters
n_params = size(particles.thetas,2);

% If a second argument isn't provided, assume all parameters are present
if nargin == 1
    present_params = 1:n_params;
end
% If a third argument isn't provided, assume no parameters are scale
% parameters
if nargin < 3
    scale_params = false(1,n_params);
else % If it was provided, just guarantee it's a logical
    scale_params = logical(scale_params);
end

% Convert the scale parameters back into their original forms (they will
% be in log form in the particles struct)
particles.thetas(:,scale_params) = exp(particles.thetas(:,scale_params));

% Initialise the parameter modes vector
param_modes = nan(1, max(present_params));
   
% Loop over all parameters
for i = 1:length(present_params{k})
    
    % Find the density at a large number of points from the minimum to
    % maximum value of this parameter. So first create a set of points to
    % evaluate at
    theta = linspace( min(particles.thetas(:,i)), max(particles.thetas(:,i)), 5000);
        
    % Create the kernel-estimated density at these points
    [p, theta] = ksdensity(particles.thetas(:,i), theta);
        
    % Find the mode and store it
    [pks, locs] = findpeaks(p, theta);
    [~, I] = sort(pks, 'descend');
    locs = locs(I);
    param_modes(present_params(i)) = locs(1);
    
end