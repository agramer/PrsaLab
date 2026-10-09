function R = circ_vmcomponents(alpha, varargin)
% CIRC_VMCOMPONENTS  Number and location of the components of a circular
% distribution with von Mises mixtures (1..Kmax components).
%
%   R = circ_vmcomponents(alpha, 'Units', 'deg')
%   R = circ_vmcomponents(alpha, 'Units', 'deg', 'Kmax', 3, 'Name', 'All neurons')
%
% For K = 1..Kmax, a mixture of K von Mises distributions is fitted by
% maximum likelihood (EM; short-EM strategy with several starts, Biernacki
% et al. 2003; kappa <= KappaMax, weights >= MinWeight). The number of
% components is selected by BIC (3K-1 parameters). In addition:
%   - excess-mass test of more than one mode (circ_modetest)
%   - likelihood-ratio test of 1 vs 2 components with parametric bootstrap
%     (circ_vmmixtest)
%   - bootstrap 95% CIs of the locations of the selected components
%     (neurons resampled, K-component model refitted, components matched to
%     the original fit by minimal total angular distance)
% and, optionally, a figure of the fitted components.
%
% Options: 'Units' ('rad'|'deg'), 'Kmax' (3), 'Bmix' (500), 'Bmode' (2000),
%          'Bboot' (500; 0 = no bootstrap, CIs NaN), 'nStarts' (10), 'KappaMax' (100), 'MinWeight'
%          (0.02), 'Plot' (true), 'Polar' (true: add a polar view of the
%          histogram and fit), 'BinWidth' (10: histogram bin width in deg,
%          bins centred on 0, BinWidth, 2*BinWidth, ...; must divide 360),
%          'Tests' (true: excess-mass and 1 vs 2 likelihood-ratio tests;
%          false: fits and bootstrap CIs only, p-values NaN), 'Verbose'
%          (true: print the results), 'Name' (''), 'OutFile' ('' = do not
%          save), 'Seed' ([])
% Output:  struct with n, BIC (1..Kmax), logL, K (selected), components
%          (table: weight, location, CI, kappa; locations in the input
%          units), pExcessMass, excessMass, pMixture1vs2, the bootstrap
%          component locations (locBoot) and mean directions (meanDirBoot),
%          and the component membership of the input angles:
%            posterior     numel(alpha) x K posterior probabilities of each
%                          component of the selected model (NaN rows for
%                          NaN inputs)
%            assignment    numel(alpha) x 1, most probable component of each
%                          angle (0 = NaN input, or max posterior below
%                          'MinPosterior', default 0 = always assigned)
%            componentIdx  1 x K cell, indices into alpha of the angles
%                          assigned to each component (numbered as in
%                          components, i.e. by decreasing weight)
%            nAssigned     number of angles assigned to each component

ip = inputParser;
ip.addParameter('Units', 'rad');   ip.addParameter('Kmax', 3);
ip.addParameter('Bmix', 500);      ip.addParameter('Bmode', 2000);   ip.addParameter('Bboot', 500);
ip.addParameter('nStarts', 10);    ip.addParameter('KappaMax', 100); ip.addParameter('MinWeight', 0.02);
ip.addParameter('Plot', true);     ip.addParameter('Name', '');      ip.addParameter('OutFile', '');
ip.addParameter('Polar', true);
ip.addParameter('Tests', true);   ip.addParameter('Verbose', true);
ip.addParameter('MinPosterior', 0, @(x) isscalar(x) && x >= 0 && x < 1);
ip.addParameter('BinWidth', 10, @(x) isscalar(x) && x > 0 && abs(360/x - round(360/x)) < 1e-9);
ip.addParameter('Seed', []);
ip.parse(varargin{:});
o = ip.Results;
if ~isempty(o.Seed), rng(o.Seed); end
deg = strcmpi(o.Units, 'deg');
a = alpha(:);  valid = ~isnan(a);  a = a(valid);   % valid: positions of the angles in alpha
if deg, a = a*pi/180; end
a = mod(a, 2*pi);  n = numel(a);

% fits for K = 1..Kmax
fits = cell(o.Kmax, 1);  ll = zeros(1, o.Kmax);
for K = 1:o.Kmax
    fits{K} = fitK(a, K, o);  ll(K) = fits{K}.ll;
end
BIC = -2*ll + (3*(1:o.Kmax) - 1) * log(n);
[~, Ksel] = min(BIC);
F = fits{Ksel};
[~, ord] = sort(F.w, 'descend');
F.w = F.w(ord);  F.mu = mod(F.mu(ord), 2*pi);  F.kappa = F.kappa(ord);

% tests
pEM = NaN;  pLR = NaN;  dEM = NaN;
if o.Tests
    [pEM, dEM] = circ_modetest(a, 'B', o.Bmode);
    pLR = circ_vmmixtest(a, 'B', o.Bmix);
end

% bootstrap CIs of component locations (and mean directions of the resamples)
locB = nan(o.Bboot, Ksel);  wB = locB;  mdB = nan(o.Bboot, 1);
P = perms(1:Ksel);
for b = 1:o.Bboot
    ab = a(randi(n, n, 1));
    mdB(b) = angle(mean(exp(1i*ab)));
    fb = fitK(ab, Ksel, o);
    best = inf;
    for q = 1:size(P, 1)                          % match components to the original fit
        c = sum(abs(angle(exp(1i*(fb.mu(P(q,:)) - F.mu)))));
        if c < best, best = c;  pq = P(q,:); end
    end
    locB(b,:) = fb.mu(pq);  wB(b,:) = fb.w(pq);
end
ci = nan(Ksel, 2);                                % NaN if Bboot = 0
for j = 1:Ksel
    if o.Bboot == 0, break; end
    off = angle(exp(1i*(locB(:,j) - F.mu(j))));
    ci(j,:) = mod(F.mu(j) + prctile(off, [2.5 97.5]), 2*pi);
end

sc = 1;  if deg, sc = 180/pi; end
R.n = n;  R.logL = ll;  R.BIC = BIC;  R.K = Ksel;
R.components = table((1:Ksel)', F.w(:), F.mu(:)*sc, ci(:,1)*sc, ci(:,2)*sc, F.kappa(:), ...
    'VariableNames', {'component', 'weight', 'location', 'ciLow', 'ciHigh', 'kappa'});
R.pExcessMass = pEM;  R.excessMass = dEM;  R.pMixture1vs2 = pLR;
R.meanDir = mod(angle(mean(exp(1i*a))), 2*pi)*sc;  R.r = abs(mean(exp(1i*a)));
R.locBoot = locB*sc;  R.weightBoot = wB;  R.meanDirBoot = mdB*sc;

% component membership of each input angle: posterior probability of each
% component of the selected model; the angle is assigned to the most
% probable component (0 = NaN input, or max posterior below MinPosterior)
lp = zeros(n, Ksel);
for j = 1:Ksel, lp(:,j) = log(F.w(j)) + vmlog(a, F.mu(j), F.kappa(j)); end
post = exp(lp - max(lp, [], 2));  post = post ./ sum(post, 2);
[pmax, lab] = max(post, [], 2);
lab(pmax < o.MinPosterior) = 0;
R.posterior = nan(numel(valid), Ksel);  R.posterior(valid, :) = post;
R.assignment = zeros(numel(valid), 1);   R.assignment(valid) = lab;
R.componentIdx = arrayfun(@(j) find(R.assignment == j), 1:Ksel, 'UniformOutput', false);
R.nAssigned = cellfun(@numel, R.componentIdx);
for K = 1:o.Kmax                                   % all fitted models, locations in input units
    [~, ok] = sort(fits{K}.w, 'descend');
    R.fits(K) = struct('K', K, 'w', fits{K}.w(ok), 'location', mod(fits{K}.mu(ok), 2*pi)*sc, 'mu', fits{K}.mu(ok),...
        'kappa', fits{K}.kappa(ok), 'logL', fits{K}.ll, 'BIC', BIC(K));
end

% report
if o.Verbose
    fprintf('%s (n = %d): mean direction %.0f, r = %.3f\n', o.Name, n, mod(R.meanDir, 360*(deg) + 2*pi*(~deg)), R.r);
    fprintf('  BIC for K = 1..%d: %s  -> %d component(s)\n', o.Kmax, mat2str(round(BIC, 1)), Ksel);
    if o.Tests
        fprintf('  excess-mass test (> 1 mode): p = %.3f;  LR test 1 vs 2 components: p = %.3f\n', pEM, pLR);
    end
    for j = 1:Ksel
        fprintf('  component %d: %.0f%% at %.0f [95%% CI %.0f-%.0f], kappa %.2f\n', j, 100*F.w(j), ...
            R.components.location(j), R.components.ciLow(j), R.components.ciHigh(j), F.kappa(j));
    end
end

if o.Plot, plotFit(a, F, ci, fits, BIC, R, o); end
end

% ==============================================================================
function F = fitK(a, K, o)
n = numel(a);
if K == 1
    C = mean(cos(a));  S = mean(sin(a));
    mu = atan2(S, C);  k = min(A1inv(hypot(C, S)), o.KappaMax);
    F = struct('w', 1, 'mu', mu, 'kappa', k, 'll', sum(vmlog(a, mu, k)));
    return
end
% starts: evenly spaced around the mean direction + random data points
m = atan2(mean(sin(a)), mean(cos(a)));
mu0 = [m + (0:K-1)*2*pi/K; a(randi(n, o.nStarts, K))];
ns = size(mu0, 1);
W = ones(ns, K)/K;  Kp = 2*ones(ns, K);  Mu = mu0;
[W, Mu, Kp, ll] = em(a, W, Mu, Kp, 50, 0, o);                % short runs
[~, keep] = sort(ll, 'descend');  keep = keep(1:min(3, ns));
[W, Mu, Kp, ll] = em(a, W(keep,:), Mu(keep,:), Kp(keep,:), 5000, 1e-7, o);
[llb, i] = max(ll);
F = struct('w', W(i,:), 'mu', Mu(i,:), 'kappa', Kp(i,:), 'll', llb);
end

function [W, Mu, Kp, ll] = em(a, W, Mu, Kp, maxIt, tol, o)
n = numel(a);  [ns, K] = size(Mu);  ca = cos(a);  sa = sin(a);
llOld = -inf(1, ns);
for it = 1:maxIt + 1
    lp = zeros(n, ns, K);
    for j = 1:K
        lp(:,:,j) = log(W(:,j))' + vmlog(a, Mu(:,j)', Kp(:,j)');
    end
    mx = max(lp, [], 3);
    lse = mx + log(sum(exp(lp - mx), 3));
    ll = sum(lse, 1);
    if it > maxIt || (tol > 0 && all(ll - llOld < tol)), break; end
    llOld = ll;
    for j = 1:K
        g = exp(lp(:,:,j) - lse);                            % n x starts
        nk = sum(g, 1)';
        W(:,j) = nk / n;
        Cj = (ca' * g)';  Sj = (sa' * g)';
        Mu(:,j) = atan2(Sj, Cj);
        Kp(:,j) = min(A1inv(hypot(Cj, Sj) ./ max(nk, eps)), o.KappaMax);
    end
    W = max(W, o.MinWeight);  W = W ./ sum(W, 2);
end
end

function lp = vmlog(x, mu, k)
lp = k .* (cos(x - mu) - 1) - log(2*pi*besseli(0, k, 1));
end

function k = A1inv(R)
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

% ==============================================================================
function polarTriangle(pax, theta, r, height, width, faceColor, edgeColor)
% POLARTRIANGLE  Filled isosceles triangle in polar axes, with its apex at
% (theta, r) and pointing toward the centre. height and width are in radial
% data units. Polar axes do not support patch/fill, so the triangle is built
% in Cartesian coordinates, filled with closely spaced lines parallel to its
% base, and each line is sampled densely and converted to polar coordinates
% (a straight line in the plane is not straight in (theta, r)).
u = [cos(theta) sin(theta)];                 % radial unit vector
v = [-u(2) u(1)];                            % tangential unit vector
apex = r * u;
b1 = (r + height) * u + width/2 * v;        % base corners (outer side)
b2 = (r + height) * u - width/2 * v;
s = linspace(0, 1, 25)';                     % samples along each line
% fill: lines across the triangle at increasing distance from the apex
nFill = 40;  pts = [];
for t = linspace(0, 1, nFill)
    p1 = apex + t*(b1 - apex);  p2 = apex + t*(b2 - apex);
    pts = [pts; p1 + s.*(p2 - p1); NaN NaN]; %#ok<AGROW>
end
[thF, rF] = cart2pol(pts(:,1), pts(:,2));
polarplot(pax, thF, rF, '-', 'Color', faceColor, 'LineWidth', 1.5);
% outline
V = [apex; b1; b2; apex];  pts = [];
for e = 1:3, pts = [pts; V(e,:) + s.*(V(e+1,:) - V(e,:))]; end %#ok<AGROW>
[thO, rO] = cart2pol(pts(:,1), pts(:,2));
polarplot(pax, thO, rO, '-', 'Color', edgeColor, 'LineWidth', 0.8);
end

% ==============================================================================
function plotFit(a, F, ci, fits, BIC, R, o)
col = [0.20 0.35 0.60];
nT = 3 + o.Polar;
f = figure('Color', 'w', 'Position', [60 60 980 + 320*o.Polar 380]);
tl = tiledlayout(f, 1, nT, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, 1, [1 2]);  hold(ax, 'on');
t = (0:0.5:360)';  tr = t*pi/180;
% bins of width bw centred on 0, bw, 2*bw, ... deg (edges -bw/2, bw/2, ...)
bw = o.BinWidth;  edges = -bw/2 : bw : 360 - bw/2;
hc = histcounts(mod(a*180/pi + bw/2, 360) - bw/2, edges, 'Normalization', 'pdf') * 360;   % density relative to uniform
% linear axis 0-360: the bin centred on 0 is shown half at 0 and half at 360
bar(ax, [0:bw:360-bw 360], [hc hc(1)], 1, 'FaceColor', col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
tot = zeros(size(tr));  comp = zeros(numel(tr), numel(F.w));
for j = 1:numel(F.w)
    comp(:,j) = F.w(j) * exp(F.kappa(j)*(cos(tr - F.mu(j)) - 1)) / besseli(0, F.kappa(j), 1);
    tot = tot + comp(:,j);
    if numel(F.w) > 1, plot(ax, t, comp(:,j), '--', 'Color', col, 'LineWidth', 1); end
end
plot(ax, t, tot, '-', 'Color', col, 'LineWidth', 2.2);
yline(ax, 1, 'k:');
for j = 1:numel(F.w)
    mu = mod(F.mu(j)*180/pi, 360);  c = mod(ci(j,:)*180/pi, 360);
    yt = interp1(t, tot, mu);
    plot(ax, mu, yt + 0.1, 'v', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 7 + 6*F.w(j));
    text(ax, mu, yt + 0.3, sprintf('%.0f%s (%.0f%%)', mu, char(176), 100*F.w(j)), 'HorizontalAlignment', 'center', 'FontSize', 9);
    yb = 0.1 + 0.1*(j-1);
    if c(2) >= c(1), plot(ax, c, [yb yb], '-', 'Color', col, 'LineWidth', 3);
    else, plot(ax, [c(1) 360], [yb yb], '-', [0 c(2)], [yb yb], '-', 'Color', col, 'LineWidth', 3); end
end
set(ax, 'XLim', [0 360], 'XTick', 0:45:360, 'YLim', [0 max(2.5, max(tot) + 0.6)], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
xlabel(ax, 'Preferred direction (deg)');  ylabel(ax, 'Density (relative to uniform)');
title(ax, {sprintf('%s (n = %d): %d component%s selected by BIC', o.Name, R.n, R.K, repmat('s', 1, R.K > 1)), ...
    sprintf('excess mass (> 1 mode) p = %.3f;  1 vs 2 components p = %.3f', R.pExcessMass, R.pMixture1vs2)}, ...
    'FontWeight', 'normal', 'FontSize', 9);
% polar plot: histogram, fitted density and components (radius = density
% relative to uniform; 0 deg to the right, counterclockwise)
if o.Polar
    pax = polaraxes(tl);  pax.Layout.Tile = 3;  hold(pax, 'on');
    polarhistogram(pax, 'BinEdges', edges*pi/180, 'BinCounts', hc, ...
        'FaceColor', col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    polarplot(pax, tr, ones(size(tr)), 'k:');
    if numel(F.w) > 1
        for j = 1:numel(F.w), polarplot(pax, tr, comp(:,j), '--', 'Color', col, 'LineWidth', 1); end
    end
    polarplot(pax, tr, tot, '-', 'Color', col, 'LineWidth', 2.2);
    rmax = 1.6 * max([hc(:); tot]);
    for j = 1:numel(F.w)
        mu = F.mu(j);  c = ci(j,:);  if c(2) < c(1), c(2) = c(2) + 2*pi; end
        rt = interp1(tr, tot, mod(mu, 2*pi));
        % triangle just outside the curve, pointing toward the centre
        th = (0.07 + 0.05*F.w(j)) * rmax;                          % size grows with weight
        polarTriangle(pax, mu, rt + 0.02*rmax, th, 0.9*th, col, 'k');
        rc = rmax * (0.98 - 0.04*(j-1));                          % CI arc (outer edge)
        polarplot(pax, linspace(c(1), c(2), 60), rc*ones(1, 60), '-', 'Color', col, 'LineWidth', 3);
        % label beyond the triangle base, anchored on the side facing the
        % centre so that it extends outward
        ha = {'right', 'center', 'left'};  va = {'top', 'middle', 'bottom'};
        text(pax, mu, rt + 0.02*rmax + th + 0.03*rmax, sprintf('%.0f%s', mod(mu*180/pi, 360), char(176)), ...
            'HorizontalAlignment', ha{2 + (cos(mu) > 0.38) - (cos(mu) < -0.38)}, ...
            'VerticalAlignment', va{2 + (sin(mu) > 0.38) - (sin(mu) < -0.38)}, 'FontSize', 8);
    end
    set(pax, 'RLim', [0 rmax], 'ThetaTick', 0:45:315, 'FontSize', 8, 'RAxisLocation', 80);
    title(pax, 'Polar view (radius: density rel. to uniform)', 'FontWeight', 'normal', 'FontSize', 9);
end

ax = nexttile(tl, nT);  hold(ax, 'on');
dB = BIC - min(BIC);
bar(ax, 1:numel(BIC), dB, 0.6, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
plot(ax, R.K, 0, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
for K = 1:numel(BIC)
    text(ax, K, dB(K) + 0.06*max([dB 1]), sprintf('%.1f', dB(K)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8);
end
ylim(ax, [0 1.15*max([dB 1])]);
set(ax, 'XTick', 1:numel(BIC), 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
xlabel(ax, 'Number of components K');  ylabel(ax, '\DeltaBIC (vs best model)');
title(ax, 'Model selection', 'FontWeight', 'normal', 'FontSize', 9);
if ~isempty(o.OutFile)
    set(findall(f, 'Type', 'axestoolbar'), 'Visible', 'off');   % keep toolbars out of the export
    exportgraphics(f, [o.OutFile '.png'], 'Resolution', 300);
    savefig(f, [o.OutFile '.fig']);
end
end
