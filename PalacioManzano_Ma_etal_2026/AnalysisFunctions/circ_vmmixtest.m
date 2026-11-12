function [pval, LR, info] = circ_vmmixtest(alpha, varargin)
% CIRC_VMMIXTEST  Parametric bootstrap likelihood-ratio test of one von Mises
% distribution against a mixture of two von Mises distributions.
%
%   [pval, LR, info] = circ_vmmixtest(alpha)
%   [pval, LR, info] = circ_vmmixtest(alpha, 'Units', 'deg', 'B', 1000)
%
% H0: the angles follow one von Mises distribution VM(mu, kappa).
% H1: they follow a two-component mixture w*VM(mu1,k1) + (1-w)*VM(mu2,k2).
% Test statistic: LR = 2*(logL_mix - logL_single), both fitted by maximum
% likelihood (the mixture by EM with several starts). Because the mixture
% model is irregular under H0, LR has no chi-square distribution; its null
% distribution is obtained by a parametric bootstrap from the fitted single
% von Mises (McLachlan 1987, Appl. Stat. 36:318-324).
%
% Interpretation: rejecting H0 means that one von Mises distribution does not
% describe the data and two components fit significantly better. That is
% evidence of bimodality only if the fitted mixture density itself has two
% modes (info.nModesMixture = 2); two overlapping components can also describe
% a single asymmetric peak. Use together with circ_modetest (nonparametric).
%
% The mixture is fitted by EM: all starts run for 50 iterations, the 3 best
% are then iterated to convergence (Biernacki et al. 2003, Comput. Stat. Data
% Anal. 41:561-575). Run time: ~6 s for n = 100 and ~90 s for n = 1000 with
% B = 100 (it scales with B).
%
% Inputs:
%   alpha        vector of angles (NaNs are removed)
% Name-value options:
%   'Units'      'rad' (default) or 'deg'
%   'B'          bootstrap replicates (default 500; 0 = fit only, pval = NaN)
%   'nStarts'    random EM starts in addition to 3 deterministic ones (default 10)
%   'KappaMax'   upper bound on component concentrations (default 100, i.e.
%                circular SD ~ 5.7 deg); prevents a component from collapsing
%                onto a few identical angles, where the likelihood is unbounded
%   'MinWeight'  lower bound on the smaller mixture weight (default 0.02)
%   'Seed'       seed for the random number generator (default: not reset)
%
% Outputs:
%   pval   bootstrap p-value, (1 + #{LR* >= LR}) / (B + 1)
%   LR     likelihood-ratio statistic
%   info   struct with the fitted parameters (in the input units), log-
%          likelihoods, AIC and BIC of both models, number of modes of the
%          fitted mixture density and the bootstrap statistics

ip = inputParser;
ip.addParameter('Units', 'rad', @(s) any(strcmpi(s, {'rad', 'deg'})));
ip.addParameter('B', 500, @(x) isscalar(x) && x >= 0);
ip.addParameter('nStarts', 10, @(x) isscalar(x) && x >= 0);
ip.addParameter('KappaMax', 100, @(x) isscalar(x) && x > 0);
ip.addParameter('MinWeight', 0.02, @(x) isscalar(x) && x > 0 && x < 0.5);
ip.addParameter('Seed', [], @(x) isempty(x) || isscalar(x));
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end

theta = alpha(:);
theta = theta(~isnan(theta));
deg = strcmpi(o.Units, 'deg');
if deg, theta = theta * pi/180; end
theta = mod(theta, 2*pi);
n = numel(theta);
if n < 5, error('circ_vmmixtest: at least 5 non-missing angles are needed.'); end

% fits to the data
[mu0, k0, ll0] = fitVM(theta, o.KappaMax);
[par1, ll1] = fitMix(theta, o);
LR = max(2 * (ll1 - ll0), 0);

% parametric bootstrap under H0
LRb = zeros(o.B, 1);
for b = 1:o.B
    xb = mod(vmrnd(mu0, k0, n), 2*pi);
    [~, ~, l0b] = fitVM(xb, o.KappaMax);
    [~, l1b] = fitMix(xb, o);
    LRb(b) = max(2 * (l1b - l0b), 0);
end
pval = NaN;  if o.B > 0, pval = (1 + sum(LRb >= LR)) / (o.B + 1); end

% fitted mixture density: number of modes
t = (0:3599)' * 2*pi/3600;
f = par1.w(1) * vmpdf(t, par1.mu(1), par1.kappa(1)) + par1.w(2) * vmpdf(t, par1.mu(2), par1.kappa(2));
pk = f > circshift(f, 1) & f >= circshift(f, -1);
nModes = sum(pk);

sc = 1;  if deg, sc = 180/pi; end
info = struct('n', n, 'B', o.B, ...
    'single', struct('mu', mod(mu0, 2*pi) * sc, 'kappa', k0, 'logL', ll0, ...
                     'AIC', -2*ll0 + 2*2, 'BIC', -2*ll0 + 2*log(n)), ...
    'mixture', struct('w', par1.w, 'mu', mod(par1.mu, 2*pi) * sc, 'kappa', par1.kappa, 'logL', ll1, ...
                      'AIC', -2*ll1 + 2*5, 'BIC', -2*ll1 + 5*log(n)), ...
    'nModesMixture', nModes, 'modesMixture', t(pk) * sc, 'LRboot', LRb);
end

% ==========================================================================
function [mu, kappa, ll] = fitVM(x, kmax)
% maximum likelihood fit of one von Mises distribution
C = mean(cos(x));  S = mean(sin(x));
mu = atan2(S, C);
kappa = min(A1inv(hypot(C, S)), kmax);
ll = sum(vmlogpdf(x, mu, kappa));
end

function [best, bestll] = fitMix(x, o)
% EM for a two-component von Mises mixture. "Short EM" strategy (Biernacki
% et al. 2003): all starts run in parallel for 50 iterations, the 3 best are
% then iterated to convergence and the best final log-likelihood is kept.
n = numel(x);
C = mean(cos(x));  S = mean(sin(x));  m = atan2(S, C);
mu = [m - pi/2, m + pi/2; m, m + pi; m - pi/4, m + pi/4];      % deterministic starts
mu = [mu; x(randi(n, o.nStarts, 2))];                          % random pairs of data points
w = 0.5 * ones(size(mu));  kappa = 2 * ones(size(mu));
[w, mu, kappa, ll] = emIter(x, w, mu, kappa, 50, 0, o);        % short runs
[~, keep] = sort(ll, 'descend');  keep = keep(1:min(3, numel(keep)));
[w, mu, kappa, ll] = emIter(x, w(keep,:), mu(keep,:), kappa(keep,:), 5000, 1e-7, o);
[bestll, i] = max(ll);
best = struct('w', w(i,:), 'mu', mu(i,:), 'kappa', kappa(i,:));
end

function [w, mu, kappa, ll] = emIter(x, w, mu, kappa, maxIt, tol, o)
% EM iterations for several starts at once (rows of w, mu, kappa); stops
% when no start improves the log-likelihood by more than tol
n = numel(x);  cx = cos(x);  sx = sin(x);
llOld = -inf(1, size(mu, 1));
for it = 1:maxIt
    lp1 = log(w(:,1))' + vmlogpdf(x, mu(:,1)', kappa(:,1)');   % n x starts
    lp2 = log(w(:,2))' + vmlogpdf(x, mu(:,2)', kappa(:,2)');
    mx = max(lp1, lp2);
    lse = mx + log(exp(lp1 - mx) + exp(lp2 - mx));
    ll = sum(lse, 1);
    if tol > 0 && all(ll - llOld < tol), return; end
    llOld = ll;
    g1 = exp(lp1 - lse);  g2 = 1 - g1;                          % responsibilities
    nk = [sum(g1, 1)' sum(g2, 1)'];
    w = min(max(nk / n, o.MinWeight), 1 - o.MinWeight);
    Cj = [(cx' * g1)' (cx' * g2)'];  Sj = [(sx' * g1)' (sx' * g2)'];
    mu = atan2(Sj, Cj);
    kappa = min(A1inv(hypot(Cj, Sj) ./ max(nk, eps)), o.KappaMax);
end
% log-likelihood of the final parameters
lp1 = log(w(:,1))' + vmlogpdf(x, mu(:,1)', kappa(:,1)');
lp2 = log(w(:,2))' + vmlogpdf(x, mu(:,2)', kappa(:,2)');
mx = max(lp1, lp2);
ll = sum(mx + log(exp(lp1 - mx) + exp(lp2 - mx)), 1);
end

function k = A1inv(R)
% inverse of A1(k) = I1(k)/I0(k), elementwise: Fisher (1993) approximation
% refined by Newton steps
R = min(max(R, 0), 1 - 1e-12);
k = 1 ./ (R.^3 - 4*R.^2 + 3*R);
i1 = R < 0.53;  k(i1) = 2*R(i1) + R(i1).^3 + 5*R(i1).^5/6;
i2 = R >= 0.53 & R < 0.85;  k(i2) = -0.4 + 1.39*R(i2) + 0.43 ./ (1 - R(i2));
for it = 1:8
    k = max(k, 1e-8);
    A = besseli(1, k, 1) ./ besseli(0, k, 1);
    step = (A - R) ./ (1 - A./k - A.^2);
    k = k - step;
    if all(abs(step) < 1e-10 * max(k, 1)), break; end
end
k = max(k, 0);
end

function lp = vmlogpdf(x, mu, kappa)
% log density of VM(mu, kappa), stable for large kappa (implicit expansion
% over rows of x and columns of mu, kappa)
lp = kappa .* (cos(x - mu) - 1) - log(2*pi*besseli(0, kappa, 1));
end

function f = vmpdf(x, mu, kappa)
f = exp(vmlogpdf(x, mu, kappa));
end

function x = vmrnd(mu, kappa, n)
% von Mises random numbers (Best & Fisher 1979)
if kappa < 1e-6, x = 2*pi*rand(n, 1); return; end
a = 1 + sqrt(1 + 4*kappa^2);
b = (a - sqrt(2*a)) / (2*kappa);
r = (1 + b^2) / (2*b);
x = zeros(n, 1);  done = false(n, 1);
while ~all(done)
    m = sum(~done);
    u1 = rand(m, 1);  u2 = rand(m, 1);  u3 = rand(m, 1);
    z = cos(pi * u1);
    f = (1 + r*z) ./ (r + z);
    c = kappa * (r - f);
    acc = (c .* (2 - c) - u2 > 0) | (log(c ./ u2) + 1 - c >= 0);
    idx = find(~done);
    x(idx(acc)) = mu + sign(u3(acc) - 0.5) .* acos(f(acc));
    done(idx(acc)) = true;
end
end
