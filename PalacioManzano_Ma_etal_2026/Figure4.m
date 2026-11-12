close all; clear; clc;

load('DataFiles\Dataset_WFsingCell.mat')
addpath(fullfile(pwd, 'PlotFunctions/'));
addpath(fullfile(pwd, 'AnalysisFunctions/'));
addpath(fullfile(pwd, 'CircStat2012a/'));
Suc = load('DataFiles\UenoContours.mat');

%% Figure 4A
figure('Color','w', 'Position', [10 500 300 500], 'Name', 'Figure 4A')
h_fit = cellfun(@(subT) logical(subT.h_fit), SessCellStats.fitParams_PeakDff);
CVs_EE = SessCellStats.PeakDffCV(SessCellStats.EE);
CVs_CH = SessCellStats.PeakDffCV(~SessCellStats.EE);

kappa = cellfun(@(subT) subT.kappa, SessCellStats.fitParams_PeakDff);
kappa_EE = kappa(~h_fit & SessCellStats.EE);
kappa_CH = kappa(~h_fit & ~SessCellStats.EE);


subplot(211)
hold on
h1 = histogram(CVs_EE, 'BinEdges', 0:0.1:2.5);
h1.FaceColor = [1, 0.75, 0.8]; 
h1.EdgeColor = 'none';
h1.Normalization = 'probability';

h2 = histogram(CVs_CH, 'BinEdges', 0:0.1:2.5);
h2.Normalization = 'probability';
h2.FaceColor = 'none';
h2.EdgeColor = [0, 0, 0];      
h2.LineWidth = 1;

plot(mean(CVs_EE, 'omitnan'), 0.21, 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
line([mean(CVs_EE, 'omitnan')-std(CVs_EE, 'omitnan') mean(CVs_EE, 'omitnan')+std(CVs_EE, 'omitnan')], [0.21 0.21], 'Color', 'r')
plot(mean(CVs_CH, 'omitnan'), 0.2, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
line([mean(CVs_CH, 'omitnan')-std(CVs_CH, 'omitnan') mean(CVs_CH, 'omitnan')+std(CVs_CH, 'omitnan')], [0.2 0.2], 'Color', 'k')

hold off
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'off', 'XLim', [-0.1 2.2])
xlabel('Coefficient of variation')
[~,p] = ttest2(CVs_CH, CVs_EE);
disp(['EE vs. CH coefficient of variation (two sample t-test): p=' num2str(p)])


subplot(212)
hold on
h1 = histogram(kappa_EE, 'BinEdges', 0:0.25:8);
h1.FaceColor = [1, 0.75, 0.8]; 
h1.EdgeColor = 'none';
h1.Normalization = 'probability';

h2 = histogram(kappa_CH, 'BinEdges', 0:0.25:8);
h2.Normalization = 'probability';
h2.FaceColor = 'none';
h2.EdgeColor = [0, 0, 0];      
h2.LineWidth = 1;

plot(median(kappa_EE, 'omitnan'), 0.22, 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
line([quantile(kappa_EE, 0.25) quantile(kappa_EE, 0.75)], [0.22 0.22], 'Color', 'r')
plot(median(kappa_CH, 'omitnan'), 0.2, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
line([quantile(kappa_CH, 0.25) quantile(kappa_CH, 0.75)], [0.2 0.2], 'Color', 'k')

hold off
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'off', 'XLim', [-0.2 5])
xlimVals = get(gca, 'XLim');
xlabel('Directional selectivity')
p = ranksum(kappa_CH, kappa_EE);
disp(['EE vs. CH coefficient of variation (rank sum test): p=' num2str(p)])



%% Figure 4B

figure('Color','w', 'Position', [210 500 400 500], 'Name', 'Figure 4B left')

h_fit = cellfun(@(subT) logical(subT.h_fit), SessCellStats.fitParams_PeakDff);
PDs = cellfun(@(subT) subT.mu, SessCellStats.fitParams_PeakDff);
PDs_EE = PDs(~h_fit & SessCellStats.EE);
PDs_CH = PDs(~h_fit & ~SessCellStats.EE);

CenterXY = vertcat(SessCellStats.CenterXY{:});
CenterXY_EE = CenterXY(~h_fit & SessCellStats.EE,:);
CenterXY_CH = CenterXY(~h_fit & ~SessCellStats.EE,:);


subplot(6,1,2:6)
hold on
AnisotropyInd = plotTopography(PDs_CH, CenterXY_CH, Suc, 0, gca);
PDs_CH_anisotropic = PDs_CH(AnisotropyInd>1);
CenterXY_CH_anisotropic = CenterXY_CH(AnisotropyInd>1, :);


binWidth = 0.1;
mlEdges = -4:binWidth:0.5;
[Nml, ~] = histcounts(CenterXY_CH_anisotropic(:,1), mlEdges);

subplot(6,1,1)
hold on
bar(mlEdges(1:end-1)+binWidth/2, (Nml/sum(Nml))/binWidth, 'BarWidth', 1, 'EdgeColor', [0 0 0], 'FaceColor', [1 0.4 0])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 0], 'XTick', [])
pd=fitdist(CenterXY_CH_anisotropic(:,1), 'kernel', 'kernel', 'normal');
y=pdf(pd,-4:0.01:0.5);
plot(-4:0.01:0.5, y, 'k')
hold off

[p,~,~,~]=dipTest(CenterXY_CH_anisotropic(:,1));
disp(['CH anisotropic zone neurons medio-lateral location unimodality (Hartigan''s dip test): p=' num2str(p)]);


figure('Color','w', 'Position', [410 500 400 500], 'Name', 'Figure 4B right')

subplot(6,1,2:6)
hold on
AnisotropyInd = plotTopography(PDs_EE, CenterXY_EE, Suc, 0, gca);
PDs_EE_anisotropic = PDs_EE(AnisotropyInd>1);
CenterXY_EE_anisotropic = CenterXY_EE(AnisotropyInd>1,:);

binWidth = 0.1;
mlEdges = -4:binWidth:0.5;
[Nml, ~] = histcounts(CenterXY_EE_anisotropic(:,1), mlEdges);

subplot(6,1,1)
hold on
bar(mlEdges(1:end-1)+binWidth/2, (Nml/sum(Nml))/binWidth, 'BarWidth', 1, 'EdgeColor', [0 0 0], 'FaceColor', [1 0.4 0])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 0], 'XTick', [])
pd=fitdist(CenterXY_EE_anisotropic(:,1), 'kernel', 'kernel', 'normal');
y=pdf(pd,-4:0.01:0.5);
plot(-4:0.01:0.5, y, 'k')
hold off

[p,~,~,~]=dipTest(CenterXY_EE_anisotropic(:,1));
disp(['EE anisotropic zone neurons medio-lateral location unimodality (Hartigan''s dip test): p=' num2str(p)]);

res_CH = in_fS1_CFA(CenterXY_CH_anisotropic, Suc);
res_EE = in_fS1_CFA(CenterXY_EE_anisotropic, Suc);

[p, ~] = chi2test_MP([res_CH.nTotal res_EE.nTotal], [res_CH.nEither res_EE.nEither]);
disp(['Anisotropic zone neurons within fS1/CFA (EE: ' num2str(res_EE.nEither) '/' num2str(res_EE.nTotal) ' vs. CH: '  num2str(res_CH.nEither) '/' num2str(res_CH.nTotal) ', chi2 test): p=' num2str(p)]);

[p, ~, ~] = circ_kuiper2(PDs_CH_anisotropic, PDs_EE_anisotropic, 'Units', 'deg');
disp(['EE vs. CH preferred direction distributions of neurons within anisotropic zones (Kuiper two-sample test): p=' num2str(p)]);

alpha = circ_ang2rad(PDs_CH_anisotropic);
[pval, ~] = circ_otest(alpha);
disp(['Hodges-Ajne test (CH neurons within anisotropic zone): p=' num2str(pval), ', N=', num2str(numel(alpha))]);

alpha = circ_ang2rad(PDs_EE_anisotropic);
[pval, ~] = circ_otest(alpha);
disp(['Hodges-Ajne test (EE neurons within anisotropic zone): p=' num2str(pval), ', N=', num2str(numel(alpha))]);

[p, ~, ~] = circ_modetest(PDs_CH_anisotropic, 'Units', 'deg');
disp(['CH preferred direction distributions of neurons within anisotropic zones (excess mass test for more than one mode): p=' num2str(p)]);

[p, ~, ~] = circ_modetest(PDs_EE_anisotropic, 'Units', 'deg');
disp(['EE preferred direction distributions of neurons within anisotropic zones (excess mass test for more than one mode): p=' num2str(p)]);


%% Figure 4C
figure('Color','w', 'Position', [610 500 420 500], 'Name', 'Figure 4C')

fprintf('Evaluating von Mises mixture models for CH group. Wait a few minutes...');
R_CH = circ_vmcomponents(PDs_CH_anisotropic, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(211)
plotPDs(PDs_CH_anisotropic, gca, 1, 0, [1 0.4 0], R_CH)

fprintf('Evaluating von Mises mixture models for EE group. Wait a few minutes...');
R_EE = circ_vmcomponents(PDs_EE_anisotropic, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(212)
plotPDs(PDs_EE_anisotropic, gca, 1, 0, [1 0.4 0], R_EE)

T = analyze_component_topography(PDs_CH_anisotropic, CenterXY_CH_anisotropic, R_CH, 'Units', 'deg', 'Plot', false, 'Verbose', false);
disp(['CH spatial segregation of the two components (vs. random, permutation test): p=' num2str(T.pSegregation)]);
disp(['CH distance between centroids of the two components (vs. random, permutation test): p=' num2str(T.pCentroid)]);


%% Figure 4D
figure('Color','w', 'Position', [10 100 800 400], 'Name', 'Figure 4D')
subplot(121)
plotTopography(PDs_CH, CenterXY_CH, Suc, 1, gca);
subplot(122)
plotTopography(PDs_EE, CenterXY_EE, Suc, 1, gca);

ind = CenterXY_CH_anisotropic(:,1)<-1.75; % Lateral neurons
CenterXY_CH_anisotropic_LAT = CenterXY_CH_anisotropic(ind,:);
PDs_CH_anisotropic_LAT = PDs_CH_anisotropic(ind,:);
ind = CenterXY_CH_anisotropic(:,1)>=-1.75; % Medial neurons
CenterXY_CH_anisotropic_MED = CenterXY_CH_anisotropic(ind,:);
PDs_CH_anisotropic_MED = PDs_CH_anisotropic(ind,:);

ind = CenterXY_EE_anisotropic(:,1)<-1.75; % Lateral neurons
CenterXY_EE_anisotropic_LAT = CenterXY_EE_anisotropic(ind,:);
PDs_EE_anisotropic_LAT = PDs_EE_anisotropic(ind,:);
ind = CenterXY_EE_anisotropic(:,1)>=-1.75; % Medial neurons
CenterXY_EE_anisotropic_MED = CenterXY_EE_anisotropic(ind,:);
PDs_EE_anisotropic_MED = PDs_EE_anisotropic(ind,:);


%% Figure 4E
figure('Color','w', 'Position', [210 100 600 600], 'Name', 'Figure 4E')

fprintf('Evaluating von Mises mixture models for CH lateral group. Wait a few minutes...');
R_CH_LAT = circ_vmcomponents(PDs_CH_anisotropic_LAT, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(221)
plotPDs(PDs_CH_anisotropic_LAT, gca, 1, 0, [0 0.6 1], R_CH_LAT)

fprintf('Evaluating von Mises mixture models for CH medial group. Wait a few minutes...');
R_CH_MED = circ_vmcomponents(PDs_CH_anisotropic_MED, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(222)
plotPDs(PDs_CH_anisotropic_MED, gca, 1, 0, [0.97 0 1], R_CH_MED)

fprintf('Evaluating von Mises mixture models for EE lateral group. Wait a few minutes...');
R_EE_LAT = circ_vmcomponents(PDs_EE_anisotropic_LAT, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(223)
plotPDs(PDs_EE_anisotropic_LAT, gca, 1, 0, [0 0.6 1], R_EE_LAT)

fprintf('Evaluating von Mises mixture models for EE medial group. Wait a few minutes...');
R_EE_MED = circ_vmcomponents(PDs_EE_anisotropic_MED, 'Units', 'deg', 'Kmax', 3, 'Plot', false, 'Verbose', false);
fprintf('Done!\n');
subplot(224)
plotPDs(PDs_EE_anisotropic_MED, gca, 1, 0, [0.97 0 1], R_EE_MED)

[p, d, info] = comp_location_test2(PDs_CH_anisotropic_LAT, PDs_CH_anisotropic_MED, 'Units', 'deg', 'Component', 'dominant', 'B', 1000);
disp(['CH fS1 cluster vs. CH CFA cluster dominant component locations (bootstrap test): p=' num2str(p)]);
[p, d, info] = comp_location_test2(PDs_EE_anisotropic_LAT, PDs_EE_anisotropic_MED, 'Units', 'deg', 'Component', 'dominant', 'B', 1000);
disp(['EE fS1 cluster vs. CH CFA cluster dominant component locations (bootstrap test): p=' num2str(p)]);

[p, d, info] = sep_interaction_test2(PDs_CH_anisotropic_LAT, PDs_CH_anisotropic_MED, PDs_EE_anisotropic_LAT, PDs_EE_anisotropic_MED, 'Units', 'deg', 'B', 1000, 'Measure', 'dominant');
disp(['EE vs. CH reduction in divergence between dominant components of CFA and fS1 clusters by d=' num2str(d) ' deg, (bootstrap test): p=' num2str(p)]);

[p, d, info] = comp_location_test2(PDs_CH_anisotropic_LAT, PDs_EE_anisotropic_LAT, 'Units', 'deg', 'Component', 'dominant', 'B', 1000);
disp(['EE vs. CH fS1 cluster dominant component locations (bootstrap test): p=' num2str(p)]);
[p, d, info] = comp_location_test2(PDs_CH_anisotropic_MED, PDs_EE_anisotropic_MED, 'Units', 'deg', 'Component', 'dominant', 'B', 1000);
disp(['EE vs. CH CFA cluster dominant component locations (bootstrap test): p=' num2str(p)]);

r  = @(a) abs(mean(exp(1i*deg2rad(a))));        % mean resultant length (PDs in degrees)
dr = @(a, b) r(b) - r(a);                       % r(EE) - r(CH) comparison
[p, d, null] = perm_test2(PDs_CH_anisotropic_LAT, PDs_EE_anisotropic_LAT, dr, 'nPerm', 1000, 'Tail', 'both');
disp(['EE vs. CH fS1 cluster concentrations (mean resultant lengths, two-sided permutation test): p=' num2str(p)]);
[p, d, null] = perm_test2(PDs_CH_anisotropic_MED, PDs_EE_anisotropic_MED, dr, 'nPerm', 1000, 'Tail', 'both');
disp(['EE vs. CH CFA cluster concentrations (mean resultant lengths, two-sided permutation test): p=' num2str(p)]);