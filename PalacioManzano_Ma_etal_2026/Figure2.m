close all; clear; clc;

load('DataFiles\Dataset_WFpopAct.mat')
addpath(fullfile(pwd, 'PlotFunctions/'));
addpath(fullfile(pwd, 'AnalysisFunctions/'));
addpath(fullfile(pwd, 'CircStat2012a/'));
Suc = load('DataFiles\UenoContours.mat');


%% Figure 2B
figure('Color','w', 'Position', [10 500 800 420], 'Name', 'Figure 2B')
subplot(121)
plotMesoActMap(ACTmaps_CH_Pop, Suc, gca)
title('Conventional housing')
subplot(122)
plotMesoActMap(ACTmaps_EE_Pop, Suc, gca)
title('Enriched environment')

%% Figure 2C-D
[p_vals, obs_stats] = cluster_permutation_test(ACTmaps_CH_Mice, ACTmaps_EE_Mice, 1999, 0.05, Suc);

%% Figure 2E
plotMapAvgAct(ACTmaps_CH_Mice, ACTmaps_EE_Mice)

%% Figure 2F-G
hFig = figure('Color','w', 'Position', [10 100 600 600], 'Name', 'Figure 2F');
mouseID_ch = 8;
plotMesoDirActMap(ACTmaps_CH_Mice_DIR, Suc, mouseID_ch, hFig);

sessID_ch = 34:39;
[PDs, ~] = getMesoMapPD(ACTmaps_CH_Sessions_DIR_unregistered, sessID_ch);
figure(hFig);
subplot(3,3,5)
plotPDs(PDs, gca, 0, 1, [0 0 0])


hFig = figure('Color','w', 'Position', [210 100 600 600], 'Name', 'Figure 2G');
mouseID_ee = 8;
plotMesoDirActMap(ACTmaps_EE_Mice_DIR, Suc, mouseID_ee, hFig);

sessID_ee = 39:43;
[PDs, ~] = getMesoMapPD(ACTmaps_EE_Sessions_DIR_unregistered, sessID_ee);
figure(hFig);
subplot(3,3,5)
plotPDs(PDs, gca, 0, 1, [0 0 0])


%% Figure 2H
sessID_ch = 1:size(ACTmaps_CH_Sessions_DIR_unregistered{1,1},3);
[PDs_CH, ~] = getMesoMapPD(ACTmaps_CH_Sessions_DIR_unregistered, sessID_ch);
figure('Color','w', 'Position', [410 100 800 400], 'Name', 'Figure 2H');
subplot(121)
plotPDs(PDs_CH, gca, 1, 0, [0 0 0])

sessID_ee = 1:size(ACTmaps_EE_Sessions_DIR_unregistered{1,1},3);
[PDs_EE, ~] = getMesoMapPD(ACTmaps_EE_Sessions_DIR_unregistered, sessID_ee);
subplot(122)
plotPDs(PDs_EE, gca, 1, 0, [0.8 0.2 0])

[p, ~, ~] = wallraff_test(circ_ang2rad(PDs_CH), circ_ang2rad(PDs_EE), 0);  
disp(['EE vs. CH PD distrubution concentration (Wallraff test): p=' num2str(p)])