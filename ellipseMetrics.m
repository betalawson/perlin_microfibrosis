function metrics = ellipseMetrics(pattern)
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
% metrics:   a (row) vector of the nine metrics relating to the three
%            ellipses


% Define the three proportions of power that form the metric ellipses
power_thresh1 = 0.1;
power_thresh2 = 0.5;
power_thresh3 = 0.8;

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

% Find the cutoff power values for the requested levels
thresh_1 = Ss(find(Sss > power_thresh1 * P_tot, 1));
thresh_2 = Ss(find(Sss > power_thresh2 * P_tot, 1));
thresh_3 = Ss(find(Sss > power_thresh3 * P_tot, 1));

% Determine binary masks that denote which locations do fall within these
% thresholds
mask_1 = zeros(size(S));
mask_2 = zeros(size(S));
mask_3 = zeros(size(S));
mask_1(S > thresh_1) = 1;
mask_2(S > thresh_2) = 1;
mask_3(S > thresh_3) = 1;

% Calculate the properties of the ellipses for these shapes
E1 = regionprops(mask_1, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
E2 = regionprops(mask_2, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
E3 = regionprops(mask_3, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');

% Store these values in one big long vector
metrics = [E1.Orientation, E1.MajorAxisLength, E1.MinorAxisLength, E2.Orientation, E2.MajorAxisLength, E2.MinorAxisLength, E3.Orientation, E3.MajorAxisLength, E3.MinorAxisLength];