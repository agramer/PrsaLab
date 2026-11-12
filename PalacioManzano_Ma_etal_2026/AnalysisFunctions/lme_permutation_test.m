function results = lme_permutation_test(dataCtrl, dataExp, varargin)
% LME_PERMUTATION_TEST  Linear mixed-effects model for two conditions with
%   arbitrary, potentially unequal numbers of repeated measurements per
%   subject, supporting both paired (within-subject) and unpaired
%   (between-subject) designs.
%
%   The random-effects structure is user-controlled:
%     - A within-subject correlation term  (1|Subject:Condition)  can be
%       added to paired designs to absorb condition-specific subject offsets.
%     - A random slope for Condition  (Condition|Subject)  can be added to
%       let the Condition effect vary across subjects.
%
%   Residuals are tested for normality via Shapiro-Wilk.  When residuals
%   are non-normal, a subject-level permutation test that preserves the
%   repeated-measures (and paired) structure is used instead.
%
% =========================================================================
% SYNTAX
%   results = lme_permutation_test(dataCtrl, dataExp)
%   results = lme_permutation_test(dataCtrl, dataExp, Name, Value, ...)
%
% =========================================================================
% INPUTS
%
%   dataCtrl  - Observations for the control condition.  One of:
%               * [NctrSub x NctrSess] numeric matrix (same number of
%                 measurements per subject; NaN = missing observation)
%               * Cell array {NctrSub x 1} where dataCtrl{i} is a numeric
%                 column vector of measurements for subject i (ragged: each
%                 subject may have a different number of observations)
%
%   dataExp   - Observations for the experimental condition.  Same formats.
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE PAIRS
%
%   'Paired'              true  -> subjects are the same individuals
%                                  measured under both conditions.
%                         false -> independent conditions (default).
%
%   'WithinSubjectCorr'   true  -> add (1|Subject:Condition) to model the
%                                  within-subject correlation across
%                                  conditions (only meaningful when
%                                  Paired = true; ignored otherwise with a
%                                  warning).
%                         'auto' -> add the term automatically when
%                                   Paired = true (default).
%                         false -> never add the term.
%
%   'RandomSlope'         true  -> add a random slope for Condition,
%                                  i.e. replace (1|Subject) with
%                                  (1 + Condition|Subject), allowing the
%                                  Condition effect to vary across subjects.
%                                  Requires Paired = true (the slope is
%                                  only identifiable when each subject
%                                  appears in both conditions).  An error
%                                  is raised if Paired = false.
%                         false -> random intercept only (default).
%
%   'Alpha'               Significance level for the Condition test (default 0.05)
%
%   'SWAlpha'             Significance level for Shapiro-Wilk  (default 0.05)
%
%   'NPermutations'       Number of permutation iterations     (default 10000)
%
%   'Verbose'             Print summary to command window      (default true)
%
% =========================================================================
% OUTPUT  -  struct with fields:
%
%   .lme                fitlme model object
%   .formula            formula string passed to fitlme
%   .lme_pvalue         p-value for the Condition fixed effect
%   .lme_tstat          t-statistic for the Condition fixed effect
%   .residuals          LME residuals (column vector)
%   .sw_stat            Shapiro-Wilk W statistic
%   .sw_pvalue          Shapiro-Wilk p-value
%   .normal_residuals   logical - true when Shapiro-Wilk does NOT reject
%   .perm_pvalue        permutation p-value  (NaN when not run)
%   .perm_observed      observed test statistic used in permutation test
%   .perm_null          [NPermutations x 1] null distribution
%   .paired             logical
%   .within_subj_corr   logical - whether (1|Subject:Condition) was included
%   .random_slope       logical - whether (Condition|Subject) was included
%   .alpha              significance level
%   .decision           plain-English summary string
%
% =========================================================================
% RANDOM-EFFECTS FORMULA LOGIC
%
%   The formula is assembled from up to three random-effects terms:
%
%     Term A  (1|Subject)               -- always included (random intercept).
%                                          Absorbed into Term B when
%                                          RandomSlope is true because
%                                          (1+Condition|Subject) already
%                                          contains a random intercept.
%
%     Term B  (1+Condition|Subject)     -- replaces Term A when
%                                          RandomSlope=true.  Models
%                                          subject-specific Condition effects.
%                                          Requires Paired=true.
%
%     Term C  (1|Subject:Condition)     -- added when WithinSubjectCorr=true
%                                          (or 'auto' with Paired=true).
%                                          Absorbs residual within-subject,
%                                          within-condition variance, capturing
%                                          the cross-condition correlation that
%                                          is not already explained by Term B.
%                                          When RandomSlope=true this term is
%                                          usually redundant (the slope already
%                                          accounts for condition differences
%                                          per subject) and is omitted with a
%                                          notice.
%
%   Resulting formulas:
%
%   Paired=F, WithinSubjectCorr=F, RandomSlope=F
%     Response ~ 1 + Condition + (1|Subject)
%
%   Paired=T, WithinSubjectCorr=F, RandomSlope=F
%     Response ~ 1 + Condition + (1|Subject)
%
%   Paired=T, WithinSubjectCorr=T, RandomSlope=F   [most common paired]
%     Response ~ 1 + Condition + (1|Subject) + (1|Subject:Condition)
%
%   Paired=T, WithinSubjectCorr=F, RandomSlope=T
%     Response ~ 1 + Condition + (1+Condition|Subject)
%
%   Paired=T, WithinSubjectCorr=T, RandomSlope=T   [WithinSubjectCorr ignored]
%     Response ~ 1 + Condition + (1+Condition|Subject)
%     NOTE: WithinSubjectCorr is automatically suppressed here because the
%     random slope already parameterises the per-subject Condition deviation;
%     adding (1|Subject:Condition) on top would create an unidentifiable model.
%
% =========================================================================
% PERMUTATION TEST
%
%   Test statistic: the LME Condition t-statistic, re-estimated on each
%   permuted dataset using the same formula as the observed model.
%   This model-based statistic is valid for ragged, unbalanced designs.
%
%   Unpaired permutation: subjects (with all their observations as a unit)
%     are shuffled across conditions, destroying between-condition structure
%     while preserving within-subject dependency.
%
%   Paired permutation: for each subject independently, the ctrl/exp labels
%     are randomly swapped with probability 0.5, preserving within-subject
%     structure while breaking the systematic condition effect under H0.
%
% =========================================================================
% REQUIREMENTS
%   Statistics and Machine Learning Toolbox  (fitlme)
%
% =========================================================================
% EXAMPLES
%
%   % --- Unpaired, ragged observations ------------------------------------
%   rng(1);
%   ctrl = arrayfun(@(~) randn(randi([2 6],1),1),     1:12, 'UniformOutput',false)';
%   expt = arrayfun(@(~) randn(randi([2 6],1),1) + 1, 1:10, 'UniformOutput',false)';
%   res  = lme_permutation_test(ctrl, expt);
%
%   % --- Paired, within-subject correlation term --------------------------
%   rng(2);  N = 15;
%   ctrl = arrayfun(@(~) randn(randi([3 6],1),1),     1:N, 'UniformOutput',false)';
%   expt = arrayfun(@(~) randn(randi([3 6],1),1) + 1, 1:N, 'UniformOutput',false)';
%   res  = lme_permutation_test(ctrl, expt, 'Paired', true);
%   % WithinSubjectCorr defaults to 'auto' -> term IS added (Paired=true)
%
%   % --- Paired, random slope, no within-subject correlation term ---------
%   res  = lme_permutation_test(ctrl, expt, 'Paired', true, ...
%              'RandomSlope', true);
%
%   % --- Paired, explicitly disable within-subject correlation term -------
%   res  = lme_permutation_test(ctrl, expt, 'Paired', true, ...
%              'WithinSubjectCorr', false);
%
%   % --- Paired, uniform matrix input, random slope -----------------------
%   rng(3);
%   ctrl = randn(12, 4);
%   expt = randn(12, 4) + 0.8;
%   res  = lme_permutation_test(ctrl, expt, 'Paired', true, ...
%              'RandomSlope', true, 'NPermutations', 5000);
%
% =========================================================================

    %% ------------------------------------------------------------------ %%
    %  1.  Parse and validate inputs
    %% ------------------------------------------------------------------ %%
    ip = inputParser();
    addRequired(ip,  'dataCtrl',            @(x) isnumeric(x) || iscell(x));
    addRequired(ip,  'dataExp',             @(x) isnumeric(x) || iscell(x));
    addParameter(ip, 'Paired',              false,  @(x) islogical(x) || ismember(x,[0 1]));
    addParameter(ip, 'WithinSubjectCorr',   'auto', @(x) (islogical(x) || ismember(x,[0 1])) || ...
                                                         (ischar(x)    && strcmpi(x,'auto')));
    addParameter(ip, 'RandomSlope',         false,  @(x) islogical(x) || ismember(x,[0 1]));
    addParameter(ip, 'Alpha',               0.05,   @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(ip, 'SWAlpha',             0.05,   @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(ip, 'NPermutations',       10000,  @(x) isscalar(x) && x >= 1);
    addParameter(ip, 'Verbose',             true,   @(x) islogical(x) || ismember(x,[0 1]));
    parse(ip, dataCtrl, dataExp, varargin{:});

    paired   = logical(ip.Results.Paired);
    wscRaw   = ip.Results.WithinSubjectCorr;   % 'auto', true, or false
    rsSlope  = logical(ip.Results.RandomSlope);
    alpha    = ip.Results.Alpha;
    swAlpha  = ip.Results.SWAlpha;
    nPerm    = ip.Results.NPermutations;
    verbose  = logical(ip.Results.Verbose);

    % --- Resolve WithinSubjectCorr ----------------------------------------
    if ischar(wscRaw)
        % 'auto': add term only when Paired=true
        useWSC = paired;
    else
        useWSC = logical(wscRaw);
    end

    % --- Cross-option validation ------------------------------------------
    if rsSlope && ~paired
        error('lme_permutation_test:randomSlopeUnpaired', ...
            ['RandomSlope=true requires Paired=true.  A random slope for ' ...
             'Condition is only identifiable when each subject contributes ' ...
             'observations to both conditions.']);
    end
    if useWSC && ~paired
        warning('lme_permutation_test:wscUnpaired', ...
            ['WithinSubjectCorr=true has no effect when Paired=false ' ...
             '(subjects appear in one condition only).  Term will be omitted.']);
        useWSC = false;
    end
    if useWSC && rsSlope
        % Random slope already models per-subject Condition deviations;
        % (1|Subject:Condition) on top creates an unidentifiable model.
        warning('lme_permutation_test:wscRedundant', ...
            ['WithinSubjectCorr=true is suppressed when RandomSlope=true.  ' ...
             'The random slope (1+Condition|Subject) already captures per-subject ' ...
             'condition offsets, making (1|Subject:Condition) unidentifiable.']);
        useWSC = false;
    end

    % --- Convert data to cell arrays -------------------------------------
    cellCtrl = to_cell(dataCtrl);
    cellExp  = to_cell(dataExp);
    NctrSub  = numel(cellCtrl);
    NexpSub  = numel(cellExp);

    if paired && NctrSub ~= NexpSub
        error('lme_permutation_test:pairedMismatch', ...
            ['Paired design requires the same number of subjects in both ' ...
             'conditions, but got %d (ctrl) vs %d (exp).'], NctrSub, NexpSub);
    end

    %% ------------------------------------------------------------------ %%
    %  2.  Build long-format table
    %% ------------------------------------------------------------------ %%
    tbl = build_table(cellCtrl, cellExp, paired);

    %% ------------------------------------------------------------------ %%
    %  3.  Assemble model formula
    %
    %   Base:          Response ~ 1 + Condition
    %   + random intercept OR random slope for Subject
    %   + optional within-subject correlation term
    %% ------------------------------------------------------------------ %%
    if rsSlope
        % (1+Condition|Subject) contains a random intercept and a random
        % slope; the two are allowed to correlate (unstructured 2x2 covariance).
        re_subject = '(1 + Condition | Subject)';
    else
        re_subject = '(1 | Subject)';
    end

    if useWSC
        re_wsc  = ' + (1 | Subject:Condition)';
    else
        re_wsc  = '';
    end

    formula = ['Response ~ 1 + Condition + ' re_subject re_wsc];

    %% ------------------------------------------------------------------ %%
    %  4.  Fit the LME
    %% ------------------------------------------------------------------ %%
    lme_fit    = fitlme(tbl, formula);
    res_vec    = residuals(lme_fit);
    feCoeffs   = lme_fit.Coefficients;
    conditionRow   = find(contains(feCoeffs.Name, 'Condition_'), 1);
    if isempty(conditionRow), conditionRow = 2; end
    lme_tstat  = feCoeffs.tStat(conditionRow);
    lme_pvalue = feCoeffs.pValue(conditionRow);

    %% ------------------------------------------------------------------ %%
    %  5.  Shapiro-Wilk normality test on LME residuals
    %% ------------------------------------------------------------------ %%
    [sw_h, sw_pvalue, sw_stat] = shapiro_wilk(res_vec, swAlpha);
    normal_residuals = ~sw_h;

    %% ------------------------------------------------------------------ %%
    %  6.  Permutation test  (only when residuals are non-normal)
    %% ------------------------------------------------------------------ %%
    perm_pvalue = NaN;
    perm_null   = NaN(nPerm, 1);
    perm_obs    = NaN;

    if ~normal_residuals
        perm_obs = lme_tstat;

        if paired
            perm_null = perm_paired(cellCtrl, cellExp, formula, nPerm);
        else
            perm_null = perm_unpaired(cellCtrl, cellExp, NctrSub, formula, nPerm);
        end

        valid = ~isnan(perm_null);
        if sum(~valid) > 0
            warning('lme_permutation_test:fitFailures', ...
                '%d / %d permutation fits failed and were excluded.', ...
                sum(~valid), nPerm);
        end
        perm_pvalue = mean(abs(perm_null(valid)) >= abs(perm_obs));
    end

    %% ------------------------------------------------------------------ %%
    %  7.  Decision string
    %% ------------------------------------------------------------------ %%
    if normal_residuals
        sig      = lme_pvalue < alpha;
        decision = sprintf( ...
            'LME Condition effect: t = %.3f, p = %.4f -> %s  (SW p = %.4f, residuals normal)', ...
            lme_tstat, lme_pvalue, sig_str(sig), sw_pvalue);
    else
        sig      = perm_pvalue < alpha;
        decision = sprintf( ...
            ['Permutation test: p = %.4f -> %s  ' ...
             '(SW p = %.4f, residuals NON-normal; LME p = %.4f for reference)'], ...
            perm_pvalue, sig_str(sig), sw_pvalue, lme_pvalue);
    end

    %% ------------------------------------------------------------------ %%
    %  8.  Pack results
    %% ------------------------------------------------------------------ %%
    results.lme              = lme_fit;
    results.formula          = formula;
    results.lme_pvalue       = lme_pvalue;
    results.lme_tstat        = lme_tstat;
    results.residuals        = res_vec;
    results.sw_stat          = sw_stat;
    results.sw_pvalue        = sw_pvalue;
    results.normal_residuals = normal_residuals;
    results.perm_pvalue      = perm_pvalue;
    results.perm_observed    = perm_obs;
    results.perm_null        = perm_null;
    results.paired           = paired;
    results.within_subj_corr = useWSC;
    results.random_slope     = rsSlope;
    results.alpha            = alpha;
    results.decision         = decision;

    %% ------------------------------------------------------------------ %%
    %  9.  Verbose summary
    %% ------------------------------------------------------------------ %%
    if verbose
        nObsCtrl = sum(cellfun(@numel, cellCtrl));
        nObsExp  = sum(cellfun(@numel, cellExp));
        fprintf('\n=== lme_permutation_test summary ===\n');
        fprintf('  Design:               %s\n', ternary(paired,'PAIRED','UNPAIRED'));
        fprintf('  Control:              %d subjects, %d total obs  (mean %.1f obs/subject)\n', ...
            NctrSub, nObsCtrl, nObsCtrl/NctrSub);
        fprintf('  Experimental:         %d subjects, %d total obs  (mean %.1f obs/subject)\n', ...
            NexpSub, nObsExp,  nObsExp/NexpSub);
        fprintf('  Within-subject corr:  %s\n', ternary(useWSC, 'YES  (1|Subject:Condition)', 'no'));
        fprintf('  Random slope:         %s\n', ternary(rsSlope,'YES  (1+Condition|Subject)','no'));
        fprintf('  Formula:              %s\n', formula);
        fprintf('\n--- LME Fixed Effects ---\n');
        disp(feCoeffs(:, {'Name','Estimate','SE','tStat','pValue'}));
        fprintf('--- Normality of Residuals (Shapiro-Wilk, alpha = %.2f) ---\n', swAlpha);
        fprintf('  W = %.4f,  p = %.4f  ->  residuals are %s\n', ...
            sw_stat, sw_pvalue, ternary(normal_residuals,'NORMAL','NON-NORMAL'));
        if ~normal_residuals
            fprintf('\n--- Permutation Test (%d permutations, alpha = %.2f) ---\n', nPerm, alpha);
            fprintf('  Test statistic (LME t):  observed = %.4f\n', perm_obs);
            fprintf('  Permutation p-value    = %.4f\n', perm_pvalue);
        end
        fprintf('\nDecision: %s\n', decision);
        fprintf('=====================================\n\n');
    end
end


%% ======================================================================= %%
%  INTERNAL HELPERS
%% ======================================================================= %%

% -------------------------------------------------------------------------
% to_cell
%   Converts a numeric matrix or already-cell input into a {NSub x 1} cell
%   array where each cell holds a column vector of observations for one
%   subject.  NaN values are treated as missing and silently dropped.
% -------------------------------------------------------------------------
function C = to_cell(data)
    if isnumeric(data)
        [nSub, ~] = size(data);
        C = cell(nSub, 1);
        for i = 1:nSub
            row  = data(i, :);
            C{i} = row(~isnan(row))';
        end
    elseif iscell(data)
        C = data(:);
        for i = 1:numel(C)
            v    = C{i};
            C{i} = v(~isnan(v(:)));
        end
    else
        error('lme_permutation_test:badInput', ...
            'Data must be a numeric matrix or a cell array of numeric vectors.');
    end
    nObs = cellfun(@numel, C);
    if any(nObs == 0)
        warning('lme_permutation_test:emptySubject', ...
            '%d subject(s) have no valid observations.', sum(nObs == 0));
    end
end

% -------------------------------------------------------------------------
% build_table
%   Assembles the long-format table required by fitlme.
%   Paired=true  -> ctrl and exp subjects share numeric Subject IDs so the
%                   LME can model within-subject correlation.
%   Paired=false -> exp subjects are offset so they are distinct individuals.
% -------------------------------------------------------------------------
function tbl = build_table(cellCtrl, cellExp, paired)
    NctrSub = numel(cellCtrl);
    NexpSub = numel(cellExp);

    ctrlIDs = (1 : NctrSub)';
    if paired
        expIDs = ctrlIDs;
    else
        expIDs = (NctrSub + 1 : NctrSub + NexpSub)';
    end

    [Y_c, Sub_c, Cond_c] = unpack(cellCtrl, ctrlIDs, 0);
    [Y_e, Sub_e, Cond_e] = unpack(cellExp,  expIDs,  1);

    Y    = [Y_c;      Y_e];
    Sub  = [Sub_c;    Sub_e];
    Cond = [Cond_c;   Cond_e];

    tbl = table( ...
        Y, ...
        categorical(Cond, [0 1], {'Ctrl','Exp'}), ...
        categorical(Sub), ...
        'VariableNames', {'Response','Condition','Subject'});
end

% -------------------------------------------------------------------------
% unpack
%   Expands a cell array of observation vectors into three aligned long-
%   format column vectors: observations, subject IDs, condition label.
% -------------------------------------------------------------------------
function [Y, SubVec, CondVec] = unpack(cellData, subIDs, condLabel)
    nObs_per_sub = cellfun(@numel, cellData);
    Y       = cell2mat(cellData(:));
    SubVec  = repelem(subIDs(:), nObs_per_sub(:));
    CondVec = repmat(condLabel, numel(Y), 1);
end

% -------------------------------------------------------------------------
% perm_unpaired
%   Shuffles the pool of subjects across conditions on each iteration.
% -------------------------------------------------------------------------
function null_dist = perm_unpaired(cellCtrl, cellExp, NctrSub, formula, nPerm)
    allSubs   = [cellCtrl(:); cellExp(:)];
    NTotal    = numel(allSubs);
    null_dist = NaN(nPerm, 1);
    for k = 1:nPerm
        idx   = randperm(NTotal);
        pCtrl = allSubs(idx(1         : NctrSub));
        pExp  = allSubs(idx(NctrSub+1 : end    ));
        tbl_p = build_table(pCtrl, pExp, false);
        null_dist(k) = fit_tstat(tbl_p, formula);
    end
end

% -------------------------------------------------------------------------
% perm_paired
%   Randomly swaps ctrl/exp labels for independent subsets of subjects.
% -------------------------------------------------------------------------
function null_dist = perm_paired(cellCtrl, cellExp, formula, nPerm)
    NSub      = numel(cellCtrl);
    null_dist = NaN(nPerm, 1);
    for k = 1:nPerm
        swapMask = logical(randi([0 1], NSub, 1));
        pCtrl    = cellCtrl;
        pExp     = cellExp;
        pCtrl(swapMask) = cellExp(swapMask);
        pExp(swapMask)  = cellCtrl(swapMask);
        tbl_p    = build_table(pCtrl, pExp, true);
        null_dist(k) = fit_tstat(tbl_p, formula);
    end
end

% -------------------------------------------------------------------------
% fit_tstat
%   Fits the LME and returns the Condition t-statistic; NaN on convergence failure.
% -------------------------------------------------------------------------
function tstat = fit_tstat(tbl, formula)
    try
        lme_p = fitlme(tbl, formula, 'CheckHessian', false);
        co    = lme_p.Coefficients;
        row   = find(contains(co.Name, 'Condition_'), 1);
        if isempty(row), row = 2; end
        tstat = co.tStat(row);
    catch
        tstat = NaN;
    end
end


%% ======================================================================= %%
%  LOCAL HELPER: Shapiro-Wilk test
%  Royston (1992) approximation, valid for 3 <= n <= 5000.
%% ======================================================================= %%
function [h, p, W] = shapiro_wilk(x, alpha)
    x = x(~isnan(x));
    x = sort(x(:));
    n = numel(x);

    if n < 3
        warning('lme_permutation_test:swSampleSize', ...
            'Shapiro-Wilk requires n >= 3.  Assuming normality.');
        h = false; p = 1; W = NaN;
        return
    end
    if n > 5000
        warning('lme_permutation_test:swSampleSize', ...
            'Shapiro-Wilk approximation may be inaccurate for n > 5000.');
    end

    % Blom's approximation for expected order statistics of N(0,1)
    m  = norminv(((1:n)' - 3/8) ./ (n + 1/4));
    m2 = m .^ 2;
    c  = m ./ sqrt(sum(m2));

    if n == 3
        a = [0.7071; 0; -0.7071];
    elseif n == 4
        a = [0.6872; 0.1677; -0.1677; -0.6872];
    elseif n == 5
        a = [0.6646; 0.2413; 0; -0.2413; -0.6646];
    else
        u   = 1 / sqrt(n);
        an  = polyval([-2.706056,  4.434685, -2.071190, -0.147981,  0.221157, c(end)  ], u);
        an1 = polyval([-3.582633,  5.682633, -1.752461, -0.293660,  0.042981, c(end-1)], u);
        phi = (sum(m2) - 2*m(end)^2 - 2*m(end-1)^2) / (1 - 2*an^2 - 2*an1^2);
        a_upper = c(end-2 : -1 : ceil(n/2)+1) ./ sqrt(phi);
        a_half  = [an; an1; a_upper(:)];
        a_lower = -flipud(a_half(1:floor(n/2)));
        if mod(n, 2) == 0
            a = [a_half; a_lower];
        else
            a = [a_half; 0; a_lower];
        end
    end

    if numel(a) > n, a = a(1:n); end
    if numel(a) < n, a(end+1:n) = 0; end

    W = (a' * x)^2 / sum((x - mean(x)).^2);
    W = min(max(W, 0), 1);

    if n <= 11
        gamma_c = polyval([0.459, -2.273], n);
        mu_c    = polyval([0.0, -1.2725, 1.0521], log(n));
        sig_c   = exp(polyval([1.0308, -0.1348], log(n)));
        zval    = (-log(gamma_c - log(1 - W)) - mu_c) / sig_c;
    else
        mu_c    = polyval([-1.2725, 1.0521], log(n));
        sig_c   = exp(polyval([1.0308, -0.1348], log(n)));
        zval    = (log(1 - W) - mu_c) / sig_c;
    end
    p = 1 - normcdf(zval);
    h = (p < alpha);
end


%% ======================================================================= %%
%  LOCAL HELPERS
%% ======================================================================= %%
function out = ternary(cond, valTrue, valFalse)
    if cond, out = valTrue; else, out = valFalse; end
end

function s = sig_str(is_sig)
    if is_sig, s = 'SIGNIFICANT'; else, s = 'not significant'; end
end