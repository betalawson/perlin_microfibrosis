function metrics = calculateMetrics(pattern)
% This function calculates a set of metrics for a provided pattern. Any
% set of functions can be called from here, that each generate their own 
% metric(s) for a pattern, and then these are combined into a single row
% vector that is visible to the remainder of the code.
%
% metrics = calculateMetrics(pattern)
%
% pattern:   the pattern for which metrics are to be calculated
%
% metrics:   a (row) vector of pattern metrics


% Calculate the nine FFT ellipse-derived metrics
ellipse_metrics = ellipseMetrics(pattern);

% Calculate the Moran's I index of spatial autocorrelation
morans_I = calculateMoransI(pattern);

% Combine these into a single vector
metrics = [ellipse_metrics, morans_I];