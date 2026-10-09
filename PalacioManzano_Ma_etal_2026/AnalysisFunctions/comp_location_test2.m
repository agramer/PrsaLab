function [p, d, info] = comp_location_test2(x1, x2, varargin)
% COMP_LOCATION_TEST2  Bootstrap test comparing the location of a von Mises
% mixture component (by default the dominant one) between two circular samples.
%
%   [p, d, info] = comp_location_test2(x1, x2, 'Units', 'deg')
%   [p, d, info] = comp_location_test2(R1, R2)          % results of circ_vmcomponents
%
% Each sample is fitted with circ_vmcomponents: von Mises mixtures with
% 1..Kmax components, the number selected by BIC, and B bootstrap resamples
% of the sample, each refitted with the selected number of components and
% matched to the original components. The compared component is the one
% with the largest weight ('dominant', default) or a given component index.
% With one component, its location is the mean direction.
%
% Statistic: d = location in sample 2 minus location in sample 1, wrapped
% to (-180, 180] deg. The bootstrap distribution of d combines the
% independent resamples of the two samples (resample b of sample 1 with
% resample b of sample 2).
% 95% CI: percentiles of the bootstrap differences.
% p-value (H0: d = 0): centred bootstrap test, the proportion of bootstrap
% differences that deviate from the observed d by at least |d|,
%   p = (1 + #{|wrap(d* - d)| >= |d|}) / (B + 1).
%
% Inputs:
%   x1, x2   angle vectors, or result structs of circ_vmcomponents (their
%            bootstrap resamples are then reused; they must have the same
%            number of bootstrap resamples and the same units)
% Name-value options (used when angles are given):
%   'Units'      'rad' (default) or 'deg'
%   'Kmax'       maximum number of components (default 3)
%   'B'          bootstrap resamples per sample (default 500)
%   'Component'  'dominant' (default) or component index (1 = largest weight)
%   'Alpha'      1 - confidence level of the CI (default 0.05)
%   'Seed'       seed for the random number generator (default: not reset)
%
% Outputs:
%   p      bootstrap p-value (two-sided)
%   d      observed difference of the component locations (input units)
%   info   struct: location1, location2, ci (of d), K1, K2 (components
%          selected), weight1, weight2, dBoot (bootstrap differences), B
%
% Example:
%   [p, d, info] = comp_location_test2(PDs_CH_fS1, PDs_CH_CFA, 'Units', 'deg');

ip = inputParser;
ip.addParameter('Units', 'rad', @(s) any(strcmpi(s, {'rad', 'deg'})));
ip.addParameter('Kmax', 3);
ip.addParameter('B', 500);
ip.addParameter('Component', 'dominant');
ip.addParameter('Alpha', 0.05);
ip.addParameter('Seed', []);
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end

deg = strcmpi(o.Units, 'deg');
fitOne = @(x) circ_vmcomponents(x, 'Units', o.Units, 'Kmax', o.Kmax, 'Bboot', o.B, ...
    'Tests', false, 'Plot', false, 'Verbose', false);
if isstruct(x1), R1 = x1; else, R1 = fitOne(x1); end
if isstruct(x2), R2 = x2; else, R2 = fitOne(x2); end
if size(R1.locBoot, 1) ~= size(R2.locBoot, 1)
    error('comp_location_test2: both results need the same number of bootstrap resamples.');
end
if isstruct(x1) && isstruct(x2)                   % infer units from the stored locations
    deg = max([R1.components.location; R2.components.location]) > 2*pi + 1e-9 || strcmpi(o.Units, 'deg');
end

if ischar(o.Component) || isstring(o.Component), j = 1; else, j = o.Component; end
if j > R1.K || j > R2.K
    error('comp_location_test2: component %d does not exist (K1 = %d, K2 = %d).', j, R1.K, R2.K);
end

sc = 1;  if deg, sc = pi/180; end                 % work in radians
L1 = R1.components.location(j) * sc;   L2 = R2.components.location(j) * sc;
b1 = R1.locBoot(:, j) * sc;             b2 = R2.locBoot(:, j) * sc;
wrap = @(x) angle(exp(1i*x));

dr = wrap(L2 - L1);
dB = wrap(b2 - b1);
off = wrap(dB - dr);                              % bootstrap deviations from the observed d
ci = dr + prctile(off, 100*[o.Alpha/2, 1 - o.Alpha/2]);
p = (1 + sum(abs(off) >= abs(dr))) / (numel(dB) + 1);

d = dr / sc;
info = struct('location1', mod(L1, 2*pi)/sc, 'location2', mod(L2, 2*pi)/sc, 'ci', ci/sc, ...
    'K1', R1.K, 'K2', R2.K, 'weight1', R1.components.weight(j), 'weight2', R2.components.weight(j), ...
    'dBoot', dB/sc, 'B', numel(dB));
end
