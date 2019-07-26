function D = ellipseDiscrepancy(metrics, target_metrics, invC)
% This function calculates a measure of discrepancy between a set of
% metrics associated with a pattern, and the metrics associated with a
% target pattern. The latter is distinguished because the eccentricity of
% the ellipses in the target pattern is used to further weight the
% distances calculated for angles. If the optional third argument is
% provided, the distance is the Mahalanobis distance, otherwise it is
% Euclidean (except with the period nature of angles taken into account).
%
% INPUTS:
%
% metrics:         a row vector (or matrix) of metrics for a single 
%                  (multiple) pattern(s)
% target_metrics:  the metrics for the target pattern
% (invC):          inverse covariance matrix for Mahalanobis distance
%
% OUTPUT:
%
% D:               discrepancies for each pattern

% First, calculate the base distances between the metrics
dM = metrics - target_metrics;

% The metrics vectors are of length (3 x <number of ellipses>)
% Because angle metrics have special handling, a loop is used to modify all
% of these. Structure is   < orientation, major_axis_length, minor_axis_length >
for k = 1:3:length(target_metrics)
    
    % Angles are periodic, so take this into account when calculating the
    % distances between angles
    dM(:,k) = angle_dists( metrics(:,k), target_metrics(k) );
    
    % Use the eccentricities of ellipses in the target pattern as scaling 
    % factors for the discrepancies in angle metrics
    dM(k)  = dM(k) * abs( log( target_metrics(k+1) / target_metrics(k+2) ) );

end

% If an (inverse) covariance matrix has been supplied, calculate the
% distance using Mahalanobis distance, otherwise use Euclidean
if nargin > 2
    D = sqrt( sum((dM * invC) .* dM, 2) );   % Uses a sum trick so as to only calculate diagonal elements of what would otherwise be two full matrix products
else
    D = vecnorm(dM',2)'; 
end


function dists = angle_dists(angles1, angle2)
% This function calculates the distance between two angles, taking into
% account the periodic nature of alignment angles (i.e. 89 degrees is two
% degrees away from -89 degrees)

% Store indices of angles where periodicity needs to be taken into account
p = abs( angles1 - angle2 ) > 90;

% Also store indices of cases where angle1 is bigger than angle2
b = angles1 > angle2;

% Now calculate distances using these indices
dists(p&b) = angle2 + 180 - angles1(p&b);
dists(p&~b) = angles1(p&~b) + 180 - angle2;
dists(~p) = angle2 - angles1(~p);