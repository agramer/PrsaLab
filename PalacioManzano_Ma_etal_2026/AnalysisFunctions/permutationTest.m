function [p, obsDiff, stats] = permutationTest(set1, set2, numMCReps)
% PERMUTATIONTEST2 Performs an exact or Monte Carlo 2-sample permutation test.
%
%   [p, obsDiff, stats] = permutationTest(set1, set2)
%   [p, obsDiff, stats] = permutationTest(set1, set2, numMCReps)
%
% WHY THIS IS APPROPRIATE FOR SMALL SAMPLE SIZES:
%   1. Distribution-Free (Non-Parametric): Parametric tests (like t-tests) 
%      rely on the Central Limit Theorem and assume normally distributed 
%      data—assumptions that cannot be validated with tiny samples (e.g., N1=5, N2=4). 
%      Permutation tests make no distributional assumptions.
%   2. Exact P-Values: For small combined N, the exact number of ways to divide 
%      the data into two groups (N! / (N1! * N2!)) is small enough to evaluate 
%      exhaustively. This calculates the exact empirical probability under 
%      the null hypothesis, rather than an asymptotic estimate.
%   3. Robust to Bounded Data: Fractions bounded in [0, 1] often exhibit skew 
%      or floor/ceiling effects that violate parametric assumptions. Permutation 
%      tests operate directly on the observed data points without transformation.
%
% Inputs:
%   set1, set2  - Numeric vectors (can be of different lengths N1 and N2)
%   numMCReps   - (Optional) Number of repetitions for Monte Carlo test
%                 if exact enumeration is too large. Default: 10,000.
%
% Outputs:
%   p           - Two-tailed p-value
%   obsDiff     - Observed difference in means (mean(set1) - mean(set2))
%   stats       - Struct containing detailed test metadata and distribution

    if nargin < 3 || isempty(numMCReps)
        numMCReps = 10000;
    end

    % Ensure column vectors
    set1 = set1(:);
    set2 = set2(:);
    
    n1 = numel(set1);
    n2 = numel(set2);
    N  = n1 + n2;
    data = [set1; set2];
    
    % Observed difference in means
    obsDiff = mean(set1) - mean(set2);
    S_total = sum(data);
    
    % Determine total possible combinations for exact test
    maxExactComb = 50000; % Threshold for exact enumeration
    try
        numComb = nchoosek(N, min(n1, n2));
    catch
        numComb = Inf; % Handles potential overflow for very large N
    end

    if numComb <= maxExactComb
        % --- EXACT PERMUTATION TEST ---
        % Small N allows computing ALL possible group partitions to yield
        % an exact, non-asymptotic p-value directly from the empirical data.
        stats.testType = 'Exact';
        
        % Generate all index combinations for set1
        idxG1 = nchoosek(1:N, n1); % (numComb x n1) matrix
        
        % Vectorized sum and mean differences across all permutations
        S1 = sum(data(idxG1), 2);
        permDiffs = (S1 / n1) - ((S_total - S1) / n2);
        
        stats.numPermutations = numComb;
    else
        % --- MONTE CARLO PERMUTATION TEST ---
        % Used when N is too large for full enumeration; randomly samples 
        % the permutation space to approximate the exact distribution.
        stats.testType = 'Monte Carlo';
        permDiffs = zeros(numMCReps, 1);
        
        for i = 1:numMCReps
            permIdx = randperm(N);
            S1 = sum(data(permIdx(1:n1)));
            permDiffs(i) = (S1 / n1) - ((S_total - S1) / n2);
        end
        
        stats.numPermutations = numMCReps;
    end

    % Two-tailed p-value calculation with floating-point tolerance
    epsTol = 1e-12;
    p = mean(abs(permDiffs) >= (abs(obsDiff) - epsTol));

    % Store output metadata
    stats.n1 = n1;
    stats.n2 = n2;
    stats.permDiffs = permDiffs;
end