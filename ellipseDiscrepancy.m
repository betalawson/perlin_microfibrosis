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

% Calculate the scaling factors for the angles based on the different
% eccentricities
%angle_scalings = abs( log(  [ target_metrics(2) / target_metrics(3), target_metrics(5) / target_metrics(6), target_metrics(8) / target_metrics(9) ] ) );
angle_scalings = ones(1,3);

% Calculate the base distances between the metrics
dM = metrics - target_metrics;

% Apply the periodic nature of angle metric
for k = [1,4,7]   % angle metrics are 1st, 4th, 7th
    dM(:,k) = angle_dists( metrics(:,k), target_metrics(k) );
end

% Scale the angle metrics
dM([1, 4, 7]) = dM([1, 4, 7]) .* angle_scalings;

% If an (inverse) covariance matrix has been supplied, calculate the
% distance using Mahalanobis distance, otherwise use Euclidean
if nargin > 2
    D = sqrt( sum((dM * invC) .* dM, 2) );   % Uses a sum trick so as to only calculate diagonal elements of what would otherwise be two full matrix products
else
    D = vecnorm(dM',2)'; 
end


function dists = angle_dists(angles1, angle2)
% This function calculates the distance between two angles, taking into
% account the periodic nature of alignment angles

% Store indices of angles where periodicity needs to be taken into account
p = abs( angles1 - angle2 ) > 90;

% Also store indices of cases where angle1 is bigger than angle2
b = angles1 > angle2;

% Now calculate distances using these indices
dists(p&b) = angle2 + 180 - angles1(p&b);
dists(p&~b) = angles1(p&~b) + 180 - angle2;
dists(~p) = angle2 - angles1(~p);