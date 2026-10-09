function [pval, stat, info] = circ_modetest(alpha, varargin)
% CIRC_MODETEST  Excess-mass test of unimodality for circular data.
%
%   [pval, stat, info] = circ_modetest(alpha)
%   [pval, stat, info] = circ_modetest(alpha, 'Units', 'deg', 'B', 1000)
%
% Tests H0: the circular distribution has one mode, against H1: more than one
% mode. Circular adaptation of the excess-mass test of Ameijeiras-Alonso,
% Crujeiras & Rodriguez-Casal (2019, TEST 28:900-919), method "ACR" of the R
% package multimode (function emcbw, mod0 = 1, unbounded support):
%   1. Test statistic: excess mass difference
%         Delta = sup_lambda [ E_2(lambda) - E_1(lambda) ],
%      E_k(lambda) = max over <= k disjoint arcs C_j of sum_j (P_n(C_j) - lambda*|C_j|),
%      P_n = empirical measure, |C| = arc length (radians). On the line this
%      equals 2*dip (the statistic multimode uses for mod0 = 1).
%   2. Critical bandwidth h_crit: smallest bandwidth for which the kernel
%      density estimate has one mode. A wrapped normal kernel is used because
%      it is the circular kernel for which the number of modes cannot increase
%      with the bandwidth (Huckemann et al. 2016, Ann. Stat. 44:1321-1345).
%   3. Calibration: B smoothed bootstrap samples, theta*_i = theta_J + h_crit*Z
%      (wrapped), and pval = mean(Delta < Delta*), as in multimode.
% The statistic is computed exactly (no grid over arcs or lambda).
%
% Inputs:
%   alpha      vector of angles (NaNs are removed)
% Name-value options:
%   'Units'    'rad' (default) or 'deg'
%   'B'        number of bootstrap replicates (default 500)
%   'GridSize' grid points for the kernel density estimate (default 2^12)
%   'Tol'      accuracy of the critical bandwidth in radians (default 1e-5)
%   'Seed'     seed for the random number generator (default: not reset)
%
% Outputs:
%   pval       bootstrap p-value (small values: evidence for > 1 mode)
%   stat       excess mass statistic Delta
%   info       struct with n, B, critical bandwidth (rad and deg), the
%              bootstrap statistics, the number of modes of the kernel density
%              estimate at the critical bandwidth and its mode location(s)
%
% Example:
%   x = [circ_vmrnd(0, 4, 100); circ_vmrnd(pi, 4, 100)];   % two modes
%   p = circ_modetest(x)

ip = inputParser;
ip.addParameter('Units', 'rad', @(s) any(strcmpi(s, {'rad', 'deg'})));
ip.addParameter('B', 500, @(x) isscalar(x) && x >= 1);
ip.addParameter('GridSize', 2^12, @(x) isscalar(x) && x >= 256);
ip.addParameter('Tol', 1e-5, @(x) isscalar(x) && x > 0);
ip.addParameter('Seed', [], @(x) isempty(x) || isscalar(x));
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end

theta = alpha(:);
theta = theta(~isnan(theta));
if strcmpi(o.Units, 'deg'), theta = theta * pi/180; end
theta = mod(theta, 2*pi);
n = numel(theta);
if n < 3, error('circ_modetest: at least 3 non-missing angles are needed.'); end

% tied angles give zero-length arcs; perturb them as multimode does
xs = sort(theta);
d = [diff(xs); xs(1) + 2*pi - xs(end)];
if any(d == 0)
    mind = min(d(d > 0));
    theta = mod(theta + (rand(n, 1) - 0.5) * mind, 2*pi);
    warning('circ_modetest:ties', 'Tied angles were perturbed by up to +/-%.3g rad.', mind/2);
end

% 1) test statistic
stat = excessMass(theta);

% 2) critical bandwidth (wrapped normal kernel)
hcrit = critBandwidth(theta, o.GridSize, o.Tol);

% 3) smoothed bootstrap from the kernel density estimate with h_crit
statB = zeros(o.B, 1);
for b = 1:o.B
    xb = mod(theta(randi(n, n, 1)) + hcrit * randn(n, 1), 2*pi);
    statB(b) = excessMass(xb);
end
pval = mean(stat < statB);

[nm, modeLoc] = nModes(theta, hcrit, o.GridSize);
info = struct('n', n, 'B', o.B, 'critBandwidthRad', hcrit, 'critBandwidthDeg', hcrit*180/pi, ...
    'nModesAtCrit', nm, 'modeAtCritRad', modeLoc, 'modeAtCritDeg', mod(modeLoc*180/pi, 360), ...
    'statBoot', statB);
end

% ==========================================================================
function D = excessMass(theta)
% Delta = sup_lambda [E_2(lambda) - E_1(lambda)] for circular data.
% E_1 is piecewise linear and convex in lambda; E_2 - E_1 is convex on every
% interval where E_1 is linear, so the supremum is attained at a breakpoint
% of E_1 (or in the limits lambda -> 0, where it is 0, and lambda -> inf,
% where it is 1/n). E_1 is computed from the shortest arc containing m points;
% E_2 is evaluated at the breakpoints by dynamic programming.
x = sort(theta(:));
n = numel(x);  w = 1/n;
g = [diff(x); x(1) + 2*pi - x(n)];          % g(i): gap from point i to point i+1 (circular)
G = [0; cumsum([g; g])];

% shortest arc containing m consecutive points, m = 1..n
L1 = zeros(n, 1);
if n <= 3000
    A = G((1:n)' + (0:n-1)) - G((1:n)');     % A(i,m) = arc from point i with m points
    L1 = min(A, [], 1)';
else
    for m = 2:n
        L1(m) = min(G((1:n)' + m - 1) - G((1:n)'));
    end
end

% breakpoints of E_1: slopes of the upper concave hull of (L1(m), m/n)
lam = hullSlopes(L1, (1:n)' * w);
if isempty(lam), D = w; return; end

E1 = max((1:n)' * w - L1 .* lam', [], 1)';   % E_1 at the breakpoints
E2 = twoArcs(g, w, lam);
D = max([E2 - E1; w]);
end

function lam = hullSlopes(L, P)
% slopes between consecutive vertices of the upper concave hull of (L, P),
% with P increasing and L nondecreasing
[L, idx] = sort(L);  P = P(idx);
keep = [diff(L) > 0; true];                  % for equal lengths keep the largest mass
L = L(keep);  P = P(keep);
h = zeros(numel(L), 1);  k = 0;
for i = 1:numel(L)
    while k >= 2 && (P(h(k)) - P(h(k-1))) * (L(i) - L(h(k))) <= (P(i) - P(h(k))) * (L(h(k)) - L(h(k-1)))
        k = k - 1;
    end
    k = k + 1;  h(k) = i;
end
h = h(1:k);
lam = diff(P(h)) ./ diff(L(h));
lam = lam(lam > 0 & isfinite(lam));
end

function E2 = twoArcs(g, w, lam)
% max over <= 2 disjoint arcs of sum(P_n(C) - lambda*|C|), for each lambda
% Case A: no arc contains the gap from the last to the first point.
% Case B: one arc wraps around; its complement is 1 or 2 "holes" (runs that
%         start and end with a gap), found as a minimum over holes.
n = numel(g);  lam = lam(:);
in1 = w * ones(size(lam));  in2 = -inf(size(lam));
best1 = in1;  best2 = in1;
for i = 2:n
    c = lam * g(i-1);
    in2 = max(in2 - c, best1) + w;           % extend arc 2, or start it after arc 1
    in1 = max(in1 - c, 0) + w;               % extend arc 1, or start it
    best1 = max(best1, in1);
    best2 = max(max(best2, in2), best1);
end
h1 = inf(size(lam));  h2 = h1;  m1 = h1;  m2 = h1;
for j = 1:n-1
    c = lam * g(j);
    h2 = min(h2 + w - c, m1 - c);            % hole value: points*w - lambda*gaps
    h1 = min(h1 + w - c, -c);
    m1 = min(m1, h1);
    m2 = min(m2, h2);
end
caseB = (1 - lam * 2*pi) - min(m1, m2);
E2 = max(best2, caseB);
end

% ==========================================================================
function h = critBandwidth(theta, M, tol)
% smallest wrapped-normal bandwidth with at most one mode (bisection, as bw.crit)
h0 = 1;
while nModes(theta, h0, M) <= 1
    h0 = h0 / 2;
    if h0 < 2*pi/M, h = h0; return; end
end
hf = 2 * h0;
while nModes(theta, hf, M) > 1
    hf = hf * 2;
    if hf > 1e3, h = hf; return; end
end
while hf - h0 >= tol
    hm = (h0 + hf) / 2;
    if nModes(theta, hm, M) <= 1, hf = hm; else, h0 = hm; end
end
h = hf;
end

function [nm, loc] = nModes(theta, h, M)
% number (and locations) of modes of the wrapped normal kernel density estimate
t = (0:M-1)' * 2*pi / M;
pos = theta / (2*pi) * M;
i0 = floor(pos);  fr = pos - i0;
cnt = accumarray(mod(i0, M) + 1, 1 - fr, [M 1]) + accumarray(mod(i0 + 1, M) + 1, fr, [M 1]);
k = [0:M/2, -M/2+1:-1]';
f = real(ifft(fft(cnt) .* exp(-0.5 * h^2 * k.^2)));
df = f - circshift(f, 1);                    % f(i) - f(i-1)
df(abs(df) < 1e-12 * max(abs(f))) = 0;
s = sign(df);
nz = find(s ~= 0);
if isempty(nz), nm = 0; loc = []; return; end
% carry the last nonzero sign over flat stretches (circularly)
sc = s;  last = s(nz(end));
for i = 1:M
    if sc(i) == 0, sc(i) = last; else, last = sc(i); end
end
peak = sc > 0 & circshift(sc, -1) < 0;       % rising into i, falling after i
nm = sum(peak);
loc = t(peak);
end
