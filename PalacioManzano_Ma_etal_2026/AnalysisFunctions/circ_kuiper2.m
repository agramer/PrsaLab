function [p, V, info] = circ_kuiper2(alpha1, alpha2, varargin)
% CIRC_KUIPER2  Kuiper two-sample test for circular data.
%
%   [p, V] = circ_kuiper2(alpha1, alpha2)
%   [p, V, info] = circ_kuiper2(alpha1, alpha2, 'Units', 'deg', 'nPerm', 100000)
%
% H0: both samples come from the same circular distribution. The test
% detects any difference (location, concentration or shape).
%
% Statistic: V = max(F1 - F2) + max(F2 - F1), where F1 and F2 are the
% empirical distribution functions of the two samples evaluated over the
% pooled sample. Unlike the Kolmogorov-Smirnov statistic, V does not depend
% on the choice of the zero direction, which makes it suitable for angles.
%
% p-value (asymptotic): Kuiper's limiting distribution,
%   p = 2 * sum_{j>=1} (4 j^2 lam^2 - 1) exp(-2 j^2 lam^2),
% with Stephens' (1970) small-sample correction
%   lam = (sqrt(ne) + 0.155 + 0.24/sqrt(ne)) * V,   ne = n1*n2/(n1+n2).
% Optionally, an exact permutation p-value is added ('nPerm' > 0; computed
% with perm_test2). The asymptotic p is accurate for moderate to large
% samples; use the permutation p for small samples or for reporting.
%
% Inputs:
%   alpha1, alpha2   vectors of angles (NaNs removed)
% Name-value options:
%   'Units'   'rad' (default) or 'deg'
%   'nPerm'   number of permutations (default 0 = no permutation test)
%   'Seed'    seed for the permutation test (default: not reset)
%
% Outputs:
%   p      asymptotic p-value (Stephens' correction)
%   V      Kuiper statistic
%   info   struct with n1, n2, ne, lambda, pAsymptotic and, if requested,
%          pPermutation and nPerm
%
% References: Kuiper NH (1960) Proc K Ned Akad Wet A 63:38-47; Stephens MA
% (1970) J R Stat Soc B 32:115-122; Press et al., Numerical Recipes, 14.3.

ip = inputParser;
ip.addParameter('Units', 'rad', @(s) any(strcmpi(s, {'rad', 'deg'})));
ip.addParameter('nPerm', 0, @(x) isscalar(x) && x >= 0);
ip.addParameter('Seed', [], @(x) isempty(x) || isscalar(x));
ip.parse(varargin{:});
o = ip.Results;

a1 = alpha1(:);  a1 = a1(~isnan(a1));
a2 = alpha2(:);  a2 = a2(~isnan(a2));
if strcmpi(o.Units, 'deg'), a1 = a1*pi/180;  a2 = a2*pi/180; end
a1 = mod(a1, 2*pi);  a2 = mod(a2, 2*pi);
n1 = numel(a1);  n2 = numel(a2);
if n1 < 2 || n2 < 2, error('circ_kuiper2: each sample needs at least 2 angles.'); end

V = kuiperStat(a1, a2);

ne = n1*n2 / (n1 + n2);
lam = (sqrt(ne) + 0.155 + 0.24/sqrt(ne)) * V;
p = kuiperTail(lam);

info = struct('n1', n1, 'n2', n2, 'ne', ne, 'lambda', lam, 'pAsymptotic', p);
if o.nPerm > 0
    info.pPermutation = perm_test2(a1, a2, @kuiperStat, 'nPerm', o.nPerm, 'Tail', 'right', 'Seed', o.Seed);
    info.nPerm = o.nPerm;
end
end

% ------------------------------------------------------------------------------
function V = kuiperStat(a1, a2)
% Kuiper statistic from one sort of the pooled sample; ties are handled by
% evaluating the distribution functions after the last of tied values
n1 = numel(a1);  n2 = numel(a2);
[x, ord] = sort([a1(:); a2(:)]);
g1 = [true(n1, 1); false(n2, 1)];  g1 = g1(ord);
F1 = cumsum(g1) / n1;  F2 = cumsum(~g1) / n2;
last = [diff(x) > 0; true];
d = F1(last) - F2(last);
V = max(d) + max(-d);
end

function p = kuiperTail(lam)
% P(Kuiper statistic > lam) under H0, asymptotic series
if lam < 0.4, p = 1; return; end              % series equals 1 to machine precision
j = (1:100)';
t = (4*j.^2*lam^2 - 1) .* exp(-2*j.^2*lam^2);
p = min(1, max(0, 2*sum(t)));
end
