function [p, W, z] = wallraff_test(alpha1, alpha2, degFlag)
% WALLRAFF_TEST  Non-parametric test for equal concentration of two
%   circular samples (Wallraff 1979).
%
%   [p, W, z] = wallraff_test(alpha1, alpha2)
%
%   Inputs  (angles in radians):
%     alpha1 - vector of angles, sample 1 (in radians)
%     alpha2 - vector of angles, sample 2 (in radians)
%
%   Outputs:
%     p  - two-tailed p-value
%     W  - Wilcoxon rank-sum statistic
%     z  - normal approximation of W
%
%   H0: the two populations have equal concentration.
%   Reject H0 if p < alpha (e.g. 0.05).

% References
%Batschelet, E (1981). Circular Statistics in Biology. chap. 6.10, p. 124
%Zar, J H (2010). Biostatistical analysis. sec. 27.7-8, p. 643 


if degFlag
    addpath('CircStat2012a\')
    alpha1 = circ_ang2rad(alpha1);
    alpha2 = circ_ang2rad(alpha2);
end

alpha1 = alpha1(:);   % ensure column vectors
alpha2 = alpha2(:);

n1 = numel(alpha1);
n2 = numel(alpha2);

%--- Step 1: compute each sample's mean direction ---
mu1 = angle(mean(exp(1i * alpha1)));
mu2 = angle(mean(exp(1i * alpha2)));


%--- Step 2: absolute angular deviation from the group mean ---
d1 = abs(angle(exp(1i * (alpha1 - mu1))));   % in [0, pi]
d2 = abs(angle(exp(1i * (alpha2 - mu2))));

%--- Step 3: Wilcoxon rank-sum test on the deviations ---
% ranksum tests H0: d1 and d2 come from the same distribution.
% More concentrated sample → smaller deviations → smaller ranks.
[p, ~, stats] = ranksum(d1, d2);
W = stats.ranksum;

%--- Normal approximation of W (for reporting) ---
n  = n1 + n2;
mu_W = n1 * (n + 1) / 2;
sig_W = sqrt(n1 * n2 * (n + 1) / 12);
z = (W - mu_W) / sig_W;

end