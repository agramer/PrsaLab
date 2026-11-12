close all; clear all; clc;

load('DataFiles\Dataset_WFsingCell.mat')
addpath(fullfile(pwd, 'PlotFunctions/'));
addpath(fullfile(pwd, 'AnalysisFunctions/'));
addpath(fullfile(pwd, 'CircStat2012a/'));
Suc = load('DataFiles\UenoContours.mat');


%% Figure 3D
figure('Color','w', 'Position', [410 500 400 460], 'Name', 'Figure 3D left')
h_fit = cellfun(@(subT) logical(subT.h_fit), SessCellStats.fitParams_PeakDff);

subplot(6,1,2:6)
hold on
plot(Suc.xySensory(:,1), Suc.xySensory(:,2), 'k--')
plot(Suc.xyMotor(:,1), Suc.xyMotor(:,2), 'k--')

CenterXY = vertcat(SessCellStats.CenterXY{:});

plot(CenterXY(~h_fit, 1), CenterXY(~h_fit, 2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.6 0.2], 'MarkerSize', 4)
plot(CenterXY(h_fit, 1), CenterXY(h_fit, 2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0 0], 'MarkerSize', 4)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 0], 'YLim', [-2 2], 'XGrid', 'on', 'YGrid', 'on')
xlabel('Medio-lateral position (mm)')
ylabel('Antero-posterior position (mm)')
hold off

subplot(6,1,1)
binWidth = 0.1;
mlEdges = -4:binWidth:0.5;
[Nml, ~] = histcounts(CenterXY(:,1), mlEdges);

bar(mlEdges(1:end-1)+binWidth/2, Nml, 'BarWidth', 1, 'EdgeColor', [0 0 0], 'FaceColor', [0.5 0.5 0.5])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 0], 'XTick', [])


angles_deg = 0:45:315;
hFig = figure('Color','w', 'Position', [610 600 300 300], 'Name', 'Figure 3D right top');
S = load('DataFiles\ExampleNeuron_DirTune.mat');
plotDirTune(S.trigData, 1, hFig, [0 0.6 0.2]);
[fitStats, r_obs, badFit] = fitDirTuneVonMises(S.trigData, 1);
disp(['Von Mises fit (directionally tuned neuron): mu=' num2str(fitStats(1)), ' deg, kappa=', num2str(fitStats(2)) ', p=' num2str(fitStats(3))]);

hFig = figure('Color','w', 'Position', [610 200 300 300], 'Name', 'Figure 3D right bottom');
S = load('DataFiles\ExampleNeuron_NonDir.mat');
plotDirTune(S.trigData, 1, hFig, [1 0 0]);



%% Figure 3E
figure('Color','w', 'Position', [10 300 420 420], 'Name', 'Figure 3E')
PDs = cellfun(@(subT) subT.mu, SessCellStats.fitParams_PeakDff);
PDs = PDs(~h_fit);
subplot(111)
plotPDs(PDs, gca, 1, 0, [0 0.6 0.2])
alpha = circ_ang2rad(PDs);
alpha = alpha(~isnan(alpha));
[pval, ~] = circ_rtest(alpha);
disp(['Rayleigh test (all directionally tuned neurons): p=' num2str(pval), ', N=', num2str(numel(alpha))]);


%% Figure 3F
figure('Color','w', 'Position', [210 300 900 400], 'Name', 'Figure 3F')
subplot(2,5,[1 2 6 7])
CenterXY_tuned = CenterXY(~h_fit,:);
AnisotropyInd = plotTopography(PDs, CenterXY_tuned, Suc, 0, gca);


subplot(2,5,[3 4 8 9])
hold on
plot(Suc.xySensory(:,1), Suc.xySensory(:,2), 'k--')
plot(Suc.xyMotor(:,1), Suc.xyMotor(:,2), 'k--')

plot(CenterXY_tuned(AnisotropyInd>1, 1), CenterXY_tuned(AnisotropyInd>1, 2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0.4 0], 'MarkerSize', 4)
plot(CenterXY_tuned(AnisotropyInd<=1, 1), CenterXY_tuned(AnisotropyInd<=1, 2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.2 0.6], 'MarkerSize', 4)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 0], 'YLim', [-2 2], 'XGrid', 'on', 'YGrid', 'on')
xlabel('Medio-lateral position (mm)')
ylabel('Antero-posterior position (mm)')
hold off

subplot(2,5,5)
plotPDs(PDs(AnisotropyInd<=1), gca, 1, 0, [0 0.2 0.6])
subplot(2,5,10)
plotPDs(PDs(AnisotropyInd>1), gca, 1, 0, [1 0.4 0])

alpha = circ_ang2rad(PDs(AnisotropyInd<=1));
alpha = alpha(~isnan(alpha));
[pval, ~] = circ_rtest(alpha);
disp(['Rayleigh test (neurons outside anisotropic zone): p=' num2str(pval), ', N=', num2str(numel(alpha))]);

alpha = circ_ang2rad(PDs(AnisotropyInd>1));
alpha = alpha(~isnan(alpha));
[pval, ~] = circ_rtest(alpha);
disp(['Rayleigh test (neurons within anisotropic zone): p=' num2str(pval), ', N=', num2str(numel(alpha))]);
