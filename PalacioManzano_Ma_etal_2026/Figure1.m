close all; clear; clc;

load('DataFiles\Dataset_Metabolism.mat')
addpath(fullfile(pwd, 'PlotFunctions/'));
addpath(fullfile(pwd, 'AnalysisFunctions/'));

%% Figure 1A
figure('Color','w', 'Position', [10 500 400 420], 'Name', 'Figure 1A')

idxEE = MetabData.EE;
idxCH = ~MetabData.EE;

data = MetabData.Weight;
subplot(121)
plotMetabData(data, idxEE, idxCH, [25 35], 'Weight (g)', gca)
[~,p] = ttest2(data(idxEE), data(idxCH));
disp(['CH vs. EE weight (two sample t-test): p=' num2str(p)])

data = 100*(MetabData.FatMass ./ MetabData.Weight);
subplot(122)
plotMetabData(data, idxEE, idxCH, [0 20], '% body fat', gca)
[~,p] = ttest2(data(idxEE), data(idxCH));
disp(['CH vs. EE % body fat (two sample t-test): p=' num2str(p)])


%% Figure 1B
figure('Color','w', 'Position', [210 500 380 420], 'Name', 'Figure 1B')

idxEE = MetabData.EE;
idxCH = ~MetabData.EE;

dataDay = MetabData.MetabRate_Day;
dataNight = MetabData.MetabRate_Night;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEE, idxCH, idxEE, idxCH, [0.3 0.7], 'Metabolic rate (kcal/h)', gca)
[~,p] = ttest2(dataDay(idxEE), dataDay(idxCH));
disp(['CH vs. EE day metabolic rate (two sample t-test): p=' num2str(p)])
[~,p] = ttest2(dataNight(idxEE), dataNight(idxCH));
disp(['CH vs. EE night metabolic rate (two sample t-test): p=' num2str(p)])

%% Figure 1C
figure('Color','w', 'Position', [310 500 380 420], 'Name', 'Figure 1C left')

idxEEday = Phenotype_Day.EE;
idxCHday = ~Phenotype_Day.EE;

idxEEnight = Phenotype_Night.EE;
idxCHnight = ~Phenotype_Night.EE;

dataDay = Phenotype_Day.foodCons;
dataNight = Phenotype_Night.foodCons;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0 4], 'Food consumption (g)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.foodCons(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.foodCons(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day food consumption (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.foodCons(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.foodCons(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night food consumption (LME model): p=' num2str(p_lme)])


figure('Color','w', 'Position', [410 500 380 420], 'Name', 'Figure 1C right')

dataDay = Phenotype_Day.waterCons;
dataNight = Phenotype_Night.waterCons;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0 5], 'Water consumption (g)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.waterCons(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.waterCons(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day water consumption (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.waterCons(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.waterCons(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night water consumption (LME model): p=' num2str(p_lme)])



%% Figure 1D
figure('Color','w', 'Position', [10 300 380 420], 'Name', 'Figure 1D left')

dataDay = Phenotype_Day.AvgRunBoutDist;
dataNight = Phenotype_Night.AvgRunBoutDist;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0 75], 'Mean running bout distance (m)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.AvgRunBoutDist(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.AvgRunBoutDist(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day mean run distance (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.AvgRunBoutDist(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.AvgRunBoutDist(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night mean run distance (LME model): p=' num2str(p_lme)])


figure('Color','w', 'Position', [210 300 380 420], 'Name', 'Figure 1D right')

dataDay = Phenotype_Day.AvgRunBoutVel;
dataNight = Phenotype_Night.AvgRunBoutVel;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0.1 0.4], 'Mean running bout velocity (m/s)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.AvgRunBoutVel(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.AvgRunBoutVel(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day mean run velocity (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.AvgRunBoutVel(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.AvgRunBoutVel(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night mean run velocity (LME model): p=' num2str(p_lme)])



%% Figure 1E
figure('Color','w', 'Position', [410 300 400 340], 'Name', 'Figure 1E left')

sessCH = 7;
stats = plotMechEfficiency(MetabEfficiencyData.WheelDist{sessCH,1}, MetabEfficiencyData.MetabRate{sessCH,1}, 0, 1);
plot(stats.xData, stats.yData, 'o', 'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','none')
hold on
xx = linspace(min(stats.xData), max(stats.xData), 100);
yy = stats.intercept + stats.slope * xx;
plot(xx, yy, 'LineWidth', 2, 'Color', 'k');
text(0, 0.85, sprintf('slope=%.3g kcal/h.m', stats.slope), 'Color', 'k');

sessEE = 19;
stats = plotMechEfficiency(MetabEfficiencyData.WheelDist{sessEE,1}, MetabEfficiencyData.MetabRate{sessEE,1}, 0, 1);
plot(stats.xData, stats.yData, 'o', 'MarkerFaceColor',[1 0.7 0.7], 'MarkerEdgeColor','none')
xx = linspace(min(stats.xData), max(stats.xData), 100);
yy = stats.intercept + stats.slope * xx;
plot(xx, yy, 'LineWidth', 2, 'Color', 'r');
text(30, 0.42, sprintf('slope=%.3g kcal/h.m', stats.slope), 'Color', 'r');

xlabel('Running distance (m)');
ylabel('Energy expenditure (kcal/h)');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XGrid', 'off', 'YGrid', 'off', 'XLim', [-10 200], 'YLim', [0.4 0.9])
hold off

figure('Color','w', 'Position', [610 300 400 340], 'Name', 'Figure 1E right')
subplot(111)
hold on
for sess = 1:size(MetabEfficiencyData,1)
    stats = plotMechEfficiency(MetabEfficiencyData.WheelDist{sess,1}, MetabEfficiencyData.MetabRate{sess,1}, 0, 1);
    xx = linspace(min(stats.xData), max(stats.xData), 100);
    yy = stats.intercept + stats.slope * xx;
    if MetabEfficiencyData.EE(sess) && ~isnan(MetabEfficiency.MetabEfficiency(sess,1))
        plot(xx, yy, 'LineWidth', 2, 'Color', 'r');
    elseif ~isnan(MetabEfficiency.MetabEfficiency(sess,1))
        plot(xx, yy, 'LineWidth', 2, 'Color', 'k');
    end
end
xlabel('Running distance (m)');
ylabel('Energy expenditure (kcal/h)');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XGrid', 'off', 'YGrid', 'off', 'XLim', [-10 350], 'YLim', [0.5 0.9])
hold off



%% Figure 1F-G
figure('Color','w', 'Position', [10 100 400 420], 'Name', 'Figure 1F-G')

idxEE = MetabEfficiency.EE;
idxCH = ~MetabEfficiency.EE;

data = MetabEfficiency.MetabEfficiency;
subplot(121)
plotMetabData(data, idxEE, idxCH, [0.6e-3 1.8e-3], '\DeltaFenergy / \DeltaFdistance (kcal/h.m)', gca)

% Reshape data for lme test
dataEE = MetabEfficiency.MetabEfficiency(idxEE);
nMice     = length(unique(MetabEfficiency.mouseID(idxEE)));
nSessions = length(dataEE) / nMice;
dataEE_Matrix = reshape(dataEE, nSessions, nMice)';

dataCH = MetabEfficiency.MetabEfficiency(idxCH);
nMice     = length(unique(MetabEfficiency.mouseID(idxCH)));
nSessions = length(dataCH) / nMice;
dataCH_Matrix = reshape(dataCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE metabolic efficiency (LME model): p=' num2str(p_lme)])

idxEE = MetabEfficiency.EE;
idxCH = ~MetabEfficiency.EE;

data = MetabEfficiency.MetabEfficiencyNorm;
subplot(122)
plotMetabData(data, idxEE, idxCH, [2.5e-5 6.5e-5], 'Normalized \DeltaFenergy / \DeltaFdistance (kcal/h.m.g)', gca)

% Reshape data for lme test
dataEE = MetabEfficiency.MetabEfficiencyNorm(idxEE);
nMice     = length(unique(MetabEfficiency.mouseID(idxEE)));
nSessions = length(dataEE) / nMice;
dataEE_Matrix = reshape(dataEE, nSessions, nMice)';

dataCH = MetabEfficiency.MetabEfficiencyNorm(idxCH);
nMice     = length(unique(MetabEfficiency.mouseID(idxCH)));
nSessions = length(dataCH) / nMice;
dataCH_Matrix = reshape(dataCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE normalized metabolic efficiency (LME model): p=' num2str(p_lme)])




%% Figure 1H
figure('Color','w', 'Position', [210 100 380 420], 'Name', 'Figure 1H')

dataDay = Phenotype_Day.fineAct;
dataNight = Phenotype_Night.fineAct;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [10 80], 'Grooming/scratching (m)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.fineAct(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.fineAct(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day grooming/scratching (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.fineAct(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.fineAct(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night mean grooming/scratching (LME model): p=' num2str(p_lme)])




%% Figure 1I
figure('Color','w', 'Position', [210 100 380 420], 'Name', 'Figure 1I')

dataDay = Phenotype_Day.AvgRearPerc;
dataNight = Phenotype_Night.AvgRearPerc;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0 80], 'Rearing (%)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.AvgRearPerc(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.AvgRearPerc(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day rearing (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.AvgRearPerc(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.AvgRearPerc(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night mean rearing (LME model): p=' num2str(p_lme)])




%% Figure 1J
figure('Color','w', 'Position', [410 100 380 420], 'Name', 'Figure 1J')

dataDay = Phenotype_Day.HabEntries;
dataNight = Phenotype_Night.HabEntries;
subplot(111)
plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, [0 50], 'Confined habitat use (# of entries)', gca)

% Reshape data for lme test
dataDayEE = Phenotype_Day.HabEntries(idxEEday);
nMice     = length(unique(Phenotype_Day.mouseID(idxEEday)));
nSessions = length(dataDayEE) / nMice;
dataEE_Matrix = reshape(dataDayEE, nSessions, nMice)';

dataDayCH = Phenotype_Day.HabEntries(idxCHday);
nMice     = length(unique(Phenotype_Day.mouseID(idxCHday)));
nSessions = length(dataDayCH) / nMice;
dataCH_Matrix = reshape(dataDayCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE day habitat entries (LME model): p=' num2str(p_lme)])

dataNightEE = Phenotype_Night.HabEntries(idxEEnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxEEnight)));
nSessions = length(dataNightEE) / nMice;
dataEE_Matrix = reshape(dataNightEE, nSessions, nMice)';

dataNightCH = Phenotype_Night.HabEntries(idxCHnight);
nMice     = length(unique(Phenotype_Night.mouseID(idxCHnight)));
nSessions = length(dataNightCH) / nMice;
dataCH_Matrix = reshape(dataNightCH, nSessions, nMice)';

results = lme_permutation_test(dataCH_Matrix, dataEE_Matrix, 'Paired', false, 'NPermutations', 1999, 'Verbose', false);
if results.normal_residuals
    p_lme = results.lme_pvalue;
else
    p_lme = results.perm_pvalue;
end
disp(['CH vs. EE night habitat entries (LME model): p=' num2str(p_lme)])