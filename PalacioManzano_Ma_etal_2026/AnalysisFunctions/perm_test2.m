function [p, obs, null] = perm_test2(x1, x2, statFun, varargin)
% PERM_TEST2  Two-sample permutation test for any statistic.
%
%   [p, obs, null] = perm_test2(x1, x2, statFun)
%   [p, obs, null] = perm_test2(x1, x2, statFun, 'nPerm', 100000, 'Tail', 'both')
%
% The group labels of the pooled observations are randomly permuted nPerm
% times; the statistic is recomputed for each permutation to build its null
% distribution (H0: both samples come from the same distribution).
%
% Inputs:
%   x1, x2    samples (vectors; NaNs are removed)
%   statFun   function handle, stat = statFun(sample1, sample2), returning a
%             scalar (e.g. @(a, b) circ_kuiper2(a, b) or a difference of
%             mean directions)
% Name-value options:
%   'nPerm'   number of permutations (default 10000)
%   'Tail'    'right' (default): large values are evidence against H0
%                                (e.g. Kuiper V, distances)
%             'both':  p based on |stat| (signed differences)
%             'left':  small values are evidence against H0
%   'Seed'    seed for the random number generator (default: not reset)
%
% Outputs:
%   p      permutation p-value, (1 + #{null as or more extreme than obs}) / (nPerm + 1)
%          (Phipson & Smyth 2010, Stat Appl Genet Mol Biol 9:39)
%   obs    observed statistic
%   null   statistics of the permuted samples (nPerm x 1)
%
% Example (difference of mean directions, angles in radians):
%   md = @(a) angle(mean(exp(1i*a)));
%   p = perm_test2(a1, a2, @(a, b) angle(exp(1i*(md(b) - md(a)))), 'Tail', 'both');

ip = inputParser;
ip.addParameter('nPerm', 10000, @(x) isscalar(x) && x >= 1);
ip.addParameter('Tail', 'right', @(s) any(strcmpi(s, {'right', 'left', 'both'})));
ip.addParameter('Seed', [], @(x) isempty(x) || isscalar(x));
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end

x1 = x1(:);  x1 = x1(~isnan(x1));
x2 = x2(:);  x2 = x2(~isnan(x2));
n1 = numel(x1);  x = [x1; x2];  N = numel(x);

obs = statFun(x1, x2);
null = zeros(o.nPerm, 1);
for k = 1:o.nPerm
    i = randperm(N);
    null(k) = statFun(x(i(1:n1)), x(i(n1+1:end)));
end

switch lower(o.Tail)
    case 'right', ext = null >= obs;
    case 'left',  ext = null <= obs;
    case 'both',  ext = abs(null) >= abs(obs);
end
p = (1 + sum(ext)) / (o.nPerm + 1);
end
