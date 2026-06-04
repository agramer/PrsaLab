function [pathLength, mid_indices] = getPathLength(xyzVals, tr)

xyzVals = cellfun(@(x) x(tr,:), xyzVals, 'UniformOutput', false);
xyzVals = cell2mat(xyzVals');
xyzVals = xyzVals';

xyzVals = xyzVals(~any(isnan(xyzVals), 2), :);

% Differences between consecutive positions
deltas = diff(xyzVals);  % size: (N-1) x 3

% Euclidean distance between each pair
segment_lengths  = sqrt(sum(deltas.^2, 2));  % size: (N-1) x 1

% Total path length
pathLength = sum(segment_lengths);

cumulative_lengths = [0; cumsum(segment_lengths)];  % N x 1

% Find index where cumulative length reaches half the path
half_length = pathLength / 2;
[~, mid_indices(2)] = min(abs(cumulative_lengths - half_length));

firstthird_length = pathLength / 3;
[~, mid_indices(1)] = min(abs(cumulative_lengths - firstthird_length));

lastthird_length = 2*pathLength / 3;
[~, mid_indices(3)] = min(abs(cumulative_lengths - lastthird_length));

