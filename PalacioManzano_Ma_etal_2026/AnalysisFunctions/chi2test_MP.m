function [pval1, pval2] = chi2test_MP(Ntot, Nlbl)

% NOTE: for a 2 by 2 test with small sample sizes use fishertest.m

% Ntot(1):   total number of samples of type 1
% Nlbl(1):   labeled number of samples of type 1
% Ntot(2):   total number of samples of type 2
% Nlbl(1):   labeled number of samples of type 2
% etc...

Nvars = numel(Ntot);
symbols = 'abcdefghijklmnopqrstuvwxyz';
% Method 1
x1=[];
x2=[];
for k=1:Nvars
    x1 = [x1; repmat(symbols(k),Ntot(k),1)];
    x2 = [x2; ones(Nlbl(k),1); 2*ones(Ntot(k)-Nlbl(k),1)];
end

[tbl,chi2stat,pval1] = crosstab(x1,x2);


% % Method 2
% p0 = (n1+n2) / (N1+N2);
% % Expected counts under H0 (null hypothesis)
% n10 = N1 * p0;
% n20 = N2 * p0;
% % Chi-square test, by hand
% observed = [n1 N1-n1 n2 N2-n2];
% expected = [n10 N1-n10 n20 N2-n20];
% chi2stat = sum((observed-expected).^2 ./ expected);
% p = 1 - chi2cdf(chi2stat,1)


% Method 3
% Pooled estimate of proportion
p0 = sum(Nlbl) / sum(Ntot);
% Expected counts under H0 (null hypothesis)
nExp = Ntot*p0;

% Chi-square test, by hand
observed = zeros(1,Nvars*2);
expected = zeros(1,Nvars*2);
for k=1:Nvars
    observed(1, 2*k-1:2*k) = [Nlbl(k) Ntot(k)-Nlbl(k)];
    expected(1, 2*k-1:2*k) = [nExp(k) Ntot(k)-nExp(k)];
end
[h,pval2,stats] = chi2gof(1:2*Nvars,'freq',observed,'expected',expected,'ctrs',1:2*Nvars,'nparams',Nvars);
