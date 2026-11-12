function T = analyze_component_topography(PDs, XY, varargin)
% ANALYZE_COMPONENT_TOPOGRAPHY  Do neurons of different von Mises mixture
% components (of their preferred directions) occupy different locations on
% the cortex, or are they intermingled?
%
%   T = analyze_component_topography(PDs, XY, 'Units', 'deg', 'Region', region)
%   T = analyze_component_topography(PDs, XY, R, ...)   % use a fitted mixture
%
% 1) Components: circ_vmcomponents (BIC-selected number of components), or
%    the already fitted mixture R (output of circ_vmcomponents for these same
%    PDs, in the same order); each neuron is assigned to its most probable
%    component (optionally only if its posterior >= 'MinPosterior';
%    unassigned neurons are excluded).
% 2) Tests, H0: component labels are unrelated to location (labels permuted
%    across neurons with perm_test2; nPerm permutations):
%    - centroid distance: distance between the mean locations of the two
%      largest components (large = spatially shifted)
%    - neighbourhood segregation: mean fraction of each neuron's k nearest
%      neighbours (on the cortex) that belong to the same component
%      (large = spatially clustered; intermingled neurons give the value
%      expected by chance)
% 3) If 'Region' is given (e.g. lateral/medial cluster of each neuron):
%    - component x region table with Fisher's exact test (2 x 2) or a
%      chi-square test with a permutation p-value (larger tables)
%    - neighbourhood segregation with labels permuted only within regions:
%      are the components clustered beyond the regional split?
% NOTE: labels are permuted across neurons, which are treated as independent.
% Neurons imaged in the same field (same mouse) are close together and may
% share preferred directions, which alone can produce spatial clustering of
% labels; a mouse-level (restricted) permutation requires mouse IDs.
%
% Inputs:  PDs  n x 1 preferred directions; XY  n x 2 locations (mm)
%          R    (optional) circ_vmcomponents result for PDs; if omitted or
%               empty, the mixture is fitted here
% Options: 'Units' ('deg'), 'Kmax' (3), 'kNN' (10), 'nPerm' (10000),
%          'MinPosterior' (0), 'Region' ([] or n x 1 labels),
%          'RegionNames' ({}), 'Name' (''), 'OutFile' (''), 'Seed' (1),
%          'Plot' (true; false: no figure, nothing saved to OutFile),
%          'Verbose' (true; false: no text output)

R = [];
if ~isempty(varargin) && (isstruct(varargin{1}) || isempty(varargin{1}))   % optional fitted mixture
    R = varargin{1};  varargin(1) = [];
end
ip = inputParser;
ip.addParameter('Units', 'deg');   ip.addParameter('Kmax', 3);
ip.addParameter('kNN', 10);        ip.addParameter('nPerm', 10000);
ip.addParameter('MinPosterior', 0);
ip.addParameter('Region', []);     ip.addParameter('RegionNames', {});
ip.addParameter('Name', '');       ip.addParameter('OutFile', '');
ip.addParameter('Seed', 1);
ip.addParameter('Plot', true);     ip.addParameter('Verbose', true);
ip.parse(varargin{:});
o = ip.Results;
o.Plot = toFlag(o.Plot, 'Plot');  o.Verbose = toFlag(o.Verbose, 'Verbose');
rng(o.Seed);

PDs = PDs(:);
if isempty(R)
    R = circ_vmcomponents(PDs, 'Units', o.Units, 'Kmax', o.Kmax, 'Bboot', 200, 'Tests', false, ...
        'Plot', false, 'Verbose', false, 'MinPosterior', o.MinPosterior, 'Seed', o.Seed);
else
    if ~isfield(R, 'posterior') || size(R.posterior, 1) ~= numel(PDs) || ~isequal(isnan(R.posterior(:,1)), isnan(PDs))
        error('analyze_component_topography: R must be the circ_vmcomponents result for these PDs (same neurons, same order).');
    end
    [pmax, lab] = max(R.posterior, [], 2);           % reassign with this function's MinPosterior
    lab(pmax < o.MinPosterior | isnan(pmax)) = 0;
    R.assignment = lab;
    R.componentIdx = arrayfun(@(j) find(lab == j), 1:R.K, 'UniformOutput', false);
    R.nAssigned = cellfun(@numel, R.componentIdx);
end
lab = R.assignment;                                  % 0 = NaN or unassigned
use = lab > 0 & all(~isnan(XY), 2);
if R.K < 2
    warning('analyze_component_topography: one component selected; nothing to compare.');
end
say(o.Verbose, '%s: %d component(s); assigned neurons per component: %s (unassigned %d)\n', o.Name, R.K, ...
    mat2str(R.nAssigned), sum(~use & ~isnan(PDs)));
for j = 1:R.K
    say(o.Verbose, '   component %d: %.0f deg (weight %.2f), n = %d, centroid ML %.2f, AP %.2f mm\n', j, ...
        R.components.location(j), R.components.weight(j), sum(lab == j & use), ...
        mean(XY(lab == j & use, 1)), mean(XY(lab == j & use, 2)));
end

idx = find(use);  L = lab(use);  P = XY(use, :);  n = numel(idx);
% k nearest neighbours among the included neurons (excluding the neuron itself)
D = sqrt((P(:,1) - P(:,1)').^2 + (P(:,2) - P(:,2)').^2);
D(1:n+1:end) = inf;
[~, nn] = sort(D, 2);  nn = nn(:, 1:o.kNN);

% statistics as functions of the index sets of component 1 and component 2
% (positions within the included neurons), for perm_test2
i1 = find(L == 1);  i2 = find(L == 2);
cent = @(a, b) norm(mean(P(a,:), 1) - mean(P(b,:), 1));
segr = @(a, b) segregation(a, b, nn, n);
[pC, dC, nullC] = perm_test2(i1, i2, cent, 'nPerm', o.nPerm, 'Tail', 'right');
[pS, sObs, nullS] = perm_test2(i1, i2, segr, 'nPerm', o.nPerm, 'Tail', 'right');
if R.K > 2
    say(o.Verbose, '   (tests use components 1 and 2; %d neurons of other components excluded)\n', sum(L > 2));
end
say(o.Verbose, '   centroid distance: %.2f mm (chance %.2f mm), permutation p = %s\n', dC, median(nullC), fmtp(pC, o.nPerm));
say(o.Verbose, '   same-component neighbours (k = %d): %.3f (chance %.3f), permutation p = %s\n', o.kNN, sObs, mean(nullS), fmtp(pS, o.nPerm));

T = struct('R', R, 'included', idx, 'labels', L, 'centroidDist', dC, 'pCentroid', pC, 'nullCentroid', nullC, ...
    'segregation', sObs, 'pSegregation', pS, 'nullSegregation', nullS, 'kNN', o.kNN);

% region analyses
if ~isempty(o.Region)
    reg = o.Region(:);  reg = reg(use);
    [ur, ~, rg] = unique(reg);
    names = o.RegionNames;  if isempty(names), names = arrayfun(@(u) sprintf('region %g', u), ur, 'UniformOutput', false); end
    tab = accumarray([L rg], 1, [R.K numel(ur)]);
    T.regionTable = tab;  T.regionNames = names;
    say(o.Verbose, '   component x region (rows: components; columns: %s):\n', strjoin(names, ', '));
    for j = 1:R.K, say(o.Verbose, '      component %d: %s\n', j, mat2str(tab(j,:))); end
    if isequal(size(tab), [2 2])
        [~, pF, st] = fishertest(tab);
        T.pFisher = pF;  T.oddsRatio = st.OddsRatio;
        say(o.Verbose, '   Fisher exact test: odds ratio %.2f, p = %.2g\n', st.OddsRatio, pF);
    end
    % segregation beyond the regional split: permute labels within regions
    ok = L <= 2;  nullW = zeros(o.nPerm, 1);
    for k = 1:o.nPerm
        Lp = L;
        for r = 1:numel(ur)
            m = find(rg == r & ok);  Lp(m) = Lp(m(randperm(numel(m))));
        end
        nullW(k) = segregation(find(Lp == 1), find(Lp == 2), nn, n); %#ok<FNDSB>
    end
    pW = (1 + sum(nullW >= sObs)) / (o.nPerm + 1);
    T.pSegregationWithinRegion = pW;  T.nullSegregationWithinRegion = nullW;
    say(o.Verbose, '   same-component neighbours, labels permuted within regions: chance %.3f, p = %s\n', mean(nullW), fmtp(pW, o.nPerm));
end

if o.Plot, plotTopo(T, XY, PDs, o); end
end

% ------------------------------------------------------------------------------
function say(verbose, varargin)
if verbose, fprintf(varargin{:}); end
end

function tf = toFlag(v, name)
% accept true/false, 1/0, 'true'/'false', 'on'/'off'
if (ischar(v) || isstring(v)) && any(strcmpi(v, {'true', 'on'})),       tf = true;
elseif (ischar(v) || isstring(v)) && any(strcmpi(v, {'false', 'off'})), tf = false;
elseif (islogical(v) || isnumeric(v)) && isscalar(v),                   tf = logical(v);
else, error('analyze_component_topography: option ''%s'' must be true or false.', name);
end
end

% ------------------------------------------------------------------------------
function s = segregation(a, b, nn, n)
% mean fraction of the k nearest neighbours of the labelled neurons that
% carry the same label (neurons in neither set count as different)
lab = zeros(n, 1);  lab(a) = 1;  lab(b) = 2;
m = lab > 0;
s = mean(mean(lab(nn(m,:)) == lab(m), 2));
end

function s = fmtp(p, N)
if p <= 1/(N+1) + eps, s = sprintf('< %.2g', 1/N); else, s = sprintf('%.4f', p); end
end

function plotTopo(T, XY, PDs, o)
R = T.R;  K = R.K;
cols = lines(max(K, 2));  cols(1:2,:) = [0.85 0.33 0.10; 0 0.45 0.74];
f = figure('Color', 'w', 'Position', [50 50 1250 440]);
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% (A) map of component membership
ax = nexttile(tl, [1 2]);  hold(ax, 'on');
if exist('UenoContours.mat', 'file')
    U = load('UenoContours.mat');
    plot(ax, U.xySensory(:,1), U.xySensory(:,2), 'k--', 'HandleVisibility', 'off');
    plot(ax, U.xyMotor(:,1), U.xyMotor(:,2), 'k--', 'HandleVisibility', 'off');
end
un = R.assignment == 0 & ~isnan(PDs);
if any(un), scatter(ax, XY(un,1), XY(un,2), 14, [0.6 0.6 0.6], 'DisplayName', 'unassigned'); end
for j = 1:K
    m = R.assignment == j;
    pj = R.posterior(m, j);
    scatter(ax, XY(m,1), XY(m,2), 10 + 30*pj, cols(j,:), 'filled', 'MarkerFaceAlpha', 0.75, ...
        'DisplayName', sprintf('component %d: %.0f%s (n = %d)', j, R.components.location(j), char(176), sum(m)));
    c = mean(XY(m,:), 1);
    plot(ax, c(1), c(2), 'p', 'MarkerSize', 16, 'MarkerFaceColor', cols(j,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end
axis(ax, 'equal');
pad = 0.3;  set(ax, 'XLim', [min(XY(:,1)) - pad, max(XY(:,1)) + pad], 'YLim', [min(XY(:,2)) - pad, max(XY(:,2)) + pad], ...
    'TickDir', 'out', 'FontSize', 9, 'Box', 'on');
xlabel(ax, 'ML (mm)');  ylabel(ax, 'AP (mm)');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, {sprintf('%s: neurons by component (size \\propto posterior; star: centroid)', o.Name)}, 'FontWeight', 'normal', 'FontSize', 9);

% (B) preferred directions by component
ax = nexttile(tl);  hold(ax, 'on');
pd = PDs;  if strcmpi(o.Units, 'rad'), pd = pd*180/pi; end
pd = mod(pd, 360);  edges = 0:15:360;
hc = zeros(numel(edges)-1, K);
for j = 1:K, hc(:,j) = histcounts(pd(R.assignment == j), edges)'; end
b = bar(ax, edges(1:end-1) + 7.5, hc, 1, 'stacked', 'EdgeColor', 'none');
for j = 1:K, b(j).FaceColor = cols(j,:); end
set(ax, 'XLim', [0 360], 'XTick', 0:90:360, 'YLim', [0 1.15*max(sum(hc, 2))], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
xlabel(ax, 'Preferred direction (deg)');  ylabel(ax, 'Neurons');
title(ax, 'Preferred directions by component', 'FontWeight', 'normal', 'FontSize', 9);

% (C) segregation test
ax = nexttile(tl);  hold(ax, 'on');
histogram(ax, T.nullSegregation, 40, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'Normalization', 'probability', ...
    'DisplayName', 'labels permuted');
if isfield(T, 'nullSegregationWithinRegion')
    histogram(ax, T.nullSegregationWithinRegion, 40, 'FaceColor', [0.3 0.6 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
        'Normalization', 'probability', 'DisplayName', 'permuted within regions');
end
xline(ax, T.segregation, 'r-', 'LineWidth', 2, 'DisplayName', 'observed');
set(ax, 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
xlabel(ax, sprintf('Same-component fraction of %d nearest neighbours', T.kNN));  ylabel(ax, 'Probability');
ttl = sprintf('Segregation: p %s', fmtp(T.pSegregation, numel(T.nullSegregation)));
if isfield(T, 'pSegregationWithinRegion')
    ttl = {ttl, sprintf('within regions: p %s', fmtp(T.pSegregationWithinRegion, numel(T.nullSegregationWithinRegion)))};
end
title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 9);
legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
if ~isempty(o.OutFile)
    set(findall(f, 'Type', 'axestoolbar'), 'Visible', 'off');
    exportgraphics(f, [o.OutFile '.png'], 'Resolution', 300);
    savefig(f, [o.OutFile '.fig']);
end
end
