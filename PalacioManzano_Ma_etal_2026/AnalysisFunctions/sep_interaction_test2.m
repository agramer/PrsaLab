function [p, d, info] = sep_interaction_test2(xA1, xA2, xB1, xB2, varargin)
% SEP_INTERACTION_TEST2  Bootstrap test of a group x condition interaction for
% circular locations: is the angular separation between two samples of
% group A different from the separation between the corresponding two
% samples of group B?
%
%   [p, d, info] = sep_interaction_test2(xA1, xA2, xB1, xB2, 'Units', 'deg')
%   [p, d, info] = sep_interaction_test2(RA1, RA2, RB1, RB2, 'Measure', 'dominant')
%
% Example: A = CH, B = EE; sample 1 = lateral cluster, sample 2 = medial
% cluster. The statistic is
%       d = |wrap(L_A2 - L_A1)| - |wrap(L_B2 - L_B1)|,
% the lateral-medial separation in CH minus that in EE, where L is the
% location of each sample:
%   'Measure', 'mean'      mean direction (default)
%   'Measure', 'dominant'  location of the dominant (largest-weight) component
%                          of the von Mises mixture selected by BIC
%                          (circ_vmcomponents)
% Separations lie in [0, 180] deg, so d lies in [-180, 180] deg; d > 0 means a
% larger separation in group A.
%
% Bootstrap: each of the four samples is resampled independently B times (or
% the resamples stored in circ_vmcomponents results are reused); resample b
% of the four samples gives the bootstrap interaction d*_b.
% 95% CI: percentiles of d*.
% p-value (H0: d = 0): centred bootstrap test,
%       p = (1 + #{|d* - d| >= |d|}) / (B + 1)
% (smallest attainable value 1/(B + 1)).
%
% Inputs:
%   xA1, xA2, xB1, xB2   angle vectors, or result structs of circ_vmcomponents
%                        (all four of the same kind; structs must have the same
%                        number of bootstrap resamples and the same units)
% Name-value options:
%   'Units'    'rad' (default) or 'deg' (for structs: units of the stored
%              locations; degrees are also detected automatically)
%   'Measure'  'mean' (default) or 'dominant'
%   'B'        bootstrap resamples per sample, when angles are given (500)
%   'Kmax'     maximum number of mixture components, 'dominant' only (3)
%   'Alpha'    1 - confidence level (0.05)
%   'Seed'     seed for the random number generator (default: not reset)
%
% Outputs:
%   p      bootstrap p-value (two-sided)
%   d      observed interaction (difference of separations; input units)
%   info   struct: sepA, sepB (separations), locations [A1 A2 B1 B2], ci,
%          dBoot (bootstrap interactions), B, measure
%
% Example:
%   [p, d, info] = sep_interaction_test2(PDs_CH_fS1, PDs_CH_CFA, PDs_EE_fS1, PDs_EE_CFA, 'Units', 'deg');

ip = inputParser;
ip.addParameter('Units', 'rad', @(s) any(strcmpi(s, {'rad', 'deg'})));
ip.addParameter('Measure', 'mean', @(s) any(strcmpi(s, {'mean', 'dominant'})));
ip.addParameter('B', 1000);
ip.addParameter('Kmax', 3);
ip.addParameter('Alpha', 0.05);
ip.addParameter('Seed', []);
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end

X = {xA1, xA2, xB1, xB2};
isS = cellfun(@isstruct, X);
if any(isS) && ~all(isS)
    error('sep_interaction_test2: give four angle vectors or four circ_vmcomponents results.');
end
deg = strcmpi(o.Units, 'deg');
dominant = strcmpi(o.Measure, 'dominant');

% locations (radians) and their bootstrap distributions, per sample
L = zeros(1, 4);  Lb = cell(1, 4);
if all(isS)
    if ~deg                                       % detect degrees from the stored locations
        deg = any(cellfun(@(R) R.meanDir > 2*pi + 1e-9 || any(R.components.location > 2*pi + 1e-9), X));
    end
    sc = 1;  if deg, sc = pi/180; end
    for k = 1:4
        R = X{k};
        if dominant, L(k) = R.components.location(1)*sc;  Lb{k} = R.locBoot(:,1)*sc;
        else,        L(k) = R.meanDir*sc;                 Lb{k} = R.meanDirBoot*sc; end
    end
else
    sc = 1;  if deg, sc = pi/180; end
    for k = 1:4
        a = X{k}(:);  a = mod(a(~isnan(a))*sc, 2*pi);
        if dominant
            R = circ_vmcomponents(a, 'Units', 'rad', 'Kmax', o.Kmax, 'Bboot', o.B, ...
                'Tests', false, 'Plot', false, 'Verbose', false);
            L(k) = R.components.location(1);  Lb{k} = R.locBoot(:,1);
        else
            n = numel(a);
            L(k) = angle(mean(exp(1i*a)));
            Lb{k} = angle(mean(exp(1i*a(randi(n, n, o.B))), 1))';   % bootstrap mean directions
        end
    end
end
nb = cellfun(@numel, Lb);
if any(nb ~= nb(1)) || nb(1) == 0
    error('sep_interaction_test2: all four samples need the same (non-zero) number of bootstrap resamples.');
end

wrap = @(x) angle(exp(1i*x));
sepA = abs(wrap(L(2) - L(1)));   sepB = abs(wrap(L(4) - L(3)));
dr = sepA - sepB;
dB = abs(wrap(Lb{2} - Lb{1})) - abs(wrap(Lb{4} - Lb{3}));
ci = prctile(dB, 100*[o.Alpha/2, 1 - o.Alpha/2]);
p = (1 + sum(abs(dB - dr) >= abs(dr))) / (numel(dB) + 1);

d = dr / sc;
info = struct('sepA', sepA/sc, 'sepB', sepB/sc, 'locations', mod(L, 2*pi)/sc, 'ci', ci/sc, ...
    'dBoot', dB/sc, 'B', numel(dB), 'measure', lower(o.Measure));
end
