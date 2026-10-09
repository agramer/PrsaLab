function [fitStats, r_obs, badFit] = fitDirTuneVonMises(trigData, cNum)

yDat = trigData.Stats{cNum,1}.PeakDff(:,1);

if min(yDat)<0
    yDat=yDat-min(yDat); % normalizes negative values to avoid NaNs
end

angles_deg = 0:45:315;

plotFlag = 0;
pThr = 0.05;    
% pThr = 0.01;    

[mu_deg, kappa, r_obs, p_value] = fitVonMises(angles_deg, yDat, 999, plotFlag);
badFit = p_value>pThr;


fitStats = [mu_deg kappa p_value];