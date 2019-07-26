function metrics = ellipseMetrics(pattern, power_thresholds)
% This function calculates a set of metrics for a provided pattern that
% quantify the relative positioning of features via the power spectrum.
% Smoothing is applied to the power spectrum, which in almost all cases
% results in ellipsoidal shapes when regions containing the top X% of power
% in the spectrum are taken. The orientation and dimensions of these
% ellipses are the metrics. Usage:
%
% metrics = ellipseMetrics(pattern)
%
% pattern:   the pattern for which metrics are to be calculated
%
% power_thresholds:   list of threshold values for %power in total spectrum
%
% metrics:   a (row) vector of the nine metrics relating to the three
%            ellipses

% Read out number of thresholds supplied
N_thresholds = length(power_thresholds);

% Calculate the power spectrum of the (mean-subtracted) pattern
P = abs(fftshift(fft2(pattern - mean(mean(pattern))))).^2;

% Apply smoothing in order to reduce noise
S = imgaussfilt(P,4);

% Calculate the total power contained
P_tot = sum(sum(S));

% Using a cumulative sum of power values sorted in descending order,
% thresholds for values that represent X% of the power can be found
Ss = sort(S(:),'descend');
Sss = cumsum(Ss);

% Loop over thresholds supplied, creating a mask for each, and then using
% this to calculate ellipse properties that become the metrics
metrics = [];
for k = 1:N_thresholds
    
    % Find the threshold power value that defines cutoff for top X% of
    % power
    pow_thresh = Ss(find(Sss > power_thresholds(k) * P_tot, 1));
    
    % Create a binary mask denoting which locations do fall within this
    % threshold
    mask = zeros(size(S));
    mask(S > pow_thresh) = 1;
    
    % Find properties of the best-fit ellipse for this mask
    E_props = regionprops(mask, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
    
    % Append these properties to the metrics vector
    metrics = [metrics, E_props.Orientation, E_props.MajorAxisLength, E_props.MinorAxisLength];
    
end