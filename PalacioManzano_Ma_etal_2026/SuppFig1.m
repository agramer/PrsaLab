close all; clear; clc;
warning off;

load('DataFiles\Dataset_WFsingCell.mat')
addpath(fullfile(pwd, 'PlotFunctions/'));
addpath(fullfile(pwd, 'AnalysisFunctions/'));
addpath(fullfile(pwd, 'CircStat2012a/'));
Suc = load('DataFiles\UenoContours.mat');


h_fit = cellfun(@(subT) logical(subT.h_fit), SessCellStats.fitParams_PeakDff);
idxEE = ~h_fit & SessCellStats.EE;
idxCH = ~h_fit & ~SessCellStats.EE;

PDs = cellfun(@(subT) subT.mu, SessCellStats.fitParams_PeakDff);
PDs_EE = PDs(~h_fit & SessCellStats.EE);
PDs_CH = PDs(~h_fit & ~SessCellStats.EE);

CenterXY = vertcat(SessCellStats.CenterXY{:});
CenterXY_EE = CenterXY(~h_fit & SessCellStats.EE,:);
CenterXY_CH = CenterXY(~h_fit & ~SessCellStats.EE,:);





mId = {'s2371';'s2372';'s2374';'s2375';'s2376'};
Nmice = size(mId,1);
figure;
AnisotropyInd = plotTopography(PDs_EE, CenterXY_EE, Suc, 0, gca); close;

Fract_in_fS1_CFA_EE = zeros(1, Nmice);
angDev_EE = zeros(1, Nmice);
angMean_EE = zeros(2, Nmice);





%% Figure S1A

figure('Color','w', 'Position', [200 100 1200 800], 'Name', 'Supp Figure 1A')
fprintf('EE mice \n');

for mouseId = 1:Nmice
    idxM = strcmp(SessCellStats.Mouse, mId(mouseId));
    Ntot = sum(idxM);
    idxM = idxEE & idxM;
    subMask = idxM(idxEE);
    Ntune = sum(subMask);

    fprintf('Mouse ID: %s \n', mId{mouseId});
    fprintf('Directionally tuned cells: N=%i/%i \n', Ntune, Ntot)

    PDs = PDs_EE(AnisotropyInd<=1 & subMask);
    subplot(4,Nmice,mouseId)
    plotPDs(PDs, gca, 1, 0, [0 0.2 0.6])
    alpha = circ_ang2rad(PDs);
    alpha = alpha(~isnan(alpha));
    [pval, ~] = circ_rtest(alpha);

    if pval < 0.001
        fprintf('Outside anisotropic zone: N=%i, p=%.3e (Rayleigh test) \n', numel(PDs), pval);
    else
        fprintf('Outside anisotropic zone: N=%i, p=%.3f (Rayleigh test) \n', numel(PDs), pval);
    end
    

    PDs = PDs_EE(AnisotropyInd>1 & subMask);
    subplot(4,Nmice,mouseId+Nmice)
    plotPDs(PDs, gca, 1, 0, [1 0.4 0])
    alpha = circ_ang2rad(PDs);
    alpha = alpha(~isnan(alpha));
    [pval, ~] = circ_rtest(alpha);
    if pval < 0.001
        fprintf('Within anisotropic zone: N=%i, p=%.3e (Rayleigh test) \n', numel(PDs), pval);
    else
        fprintf('Within anisotropic zone: N=%i, p=%.3f (Rayleigh test) \n', numel(PDs), pval);
    end
    angDev_EE(1,mouseId) = concent_param(alpha, 0);
    
    ind = AnisotropyInd>1 & subMask;
    inFlag_fS1 = inpolygon(CenterXY_EE(ind,1), CenterXY_EE(ind,2), Suc.xySensory(:,1), Suc.xySensory(:,2));
    inFlag_CFA = inpolygon(CenterXY_EE(ind,1), CenterXY_EE(ind,2), Suc.xyMotor(:,1), Suc.xyMotor(:,2));
    Nin = sum(inFlag_fS1 | inFlag_CFA);
    fprintf('Within fS1/CFA: N=%i/%i \n', Nin, numel(PDs))
    Fract_in_fS1_CFA_EE(1,mouseId) = Nin/numel(PDs);
    
    PDs = PDs_EE(AnisotropyInd>1 & subMask & CenterXY_EE(:,1)<-1.75);
    subplot(4,Nmice,mouseId+2*Nmice)
    plotPDs(PDs, gca, 1, 0, [0 0.6 1])
    fprintf('Lateral neurons: N=%i \n', numel(PDs))
    alpha1 = circ_ang2rad(PDs);
    alpha1 = alpha1(~isnan(alpha1));

    PDs = PDs_EE(AnisotropyInd>1 & subMask & CenterXY_EE(:,1)>=-1.75);
    subplot(4,Nmice,mouseId+3*Nmice)
    plotPDs(PDs, gca, 1, 0, [0.97 0 1])
    fprintf('Medial neurons: N=%i \n', numel(PDs))
    alpha2 = circ_ang2rad(PDs);
    alpha2 = alpha2(~isnan(alpha2));

    [pval, ~] = circ_wwtest(alpha1, alpha2);
    fprintf('Wason-Williams test for equal means (Lateral vs. Medial): p=%.2f', pval);

    angMean_EE(:,mouseId) = [circ_mean(alpha1);circ_mean(alpha2)];


    fprintf('\n\n\n');

end




mId = {'b2362';'b2363';'b2488';'b2561'};
Nmice = size(mId,1);
figure;
AnisotropyInd = plotTopography(PDs_CH, CenterXY_CH, Suc, 0, gca); close;

Fract_in_fS1_CFA_CH = zeros(1, Nmice);
angDev_CH = zeros(1, Nmice);
angMean_CH = zeros(2, Nmice);

%% Figure S1B

figure('Color','w', 'Position', [500 10 980 800], 'Name', 'Supp Figure 1B')
fprintf('CH mice \n');

for mouseId = 1:Nmice
    idxM = strcmp(SessCellStats.Mouse, mId(mouseId));
    Ntot = sum(idxM);
    idxM = idxCH & idxM;
    subMask = idxM(idxCH);
    Ntune = sum(subMask);

    fprintf('Mouse ID: %s \n', mId{mouseId});
    fprintf('Directionally tuned cells: N=%i/%i \n', Ntune, Ntot)

    PDs = PDs_CH(AnisotropyInd<=1 & subMask);
    subplot(4,Nmice,mouseId)
    plotPDs(PDs, gca, 1, 0, [0 0.2 0.6])
    alpha = circ_ang2rad(PDs);
    alpha = alpha(~isnan(alpha));
    [pval, ~] = circ_rtest(alpha);

    if pval < 0.001
        fprintf('Outside anisotropic zone: N=%i, p=%.3e (Rayleigh test) \n', numel(PDs), pval);
    else
        fprintf('Outside anisotropic zone: N=%i, p=%.3f (Rayleigh test) \n', numel(PDs), pval);
    end
    
    
    PDs = PDs_CH(AnisotropyInd>1 & subMask);
    subplot(4,Nmice,mouseId+Nmice)
    plotPDs(PDs, gca, 1, 0, [1 0.4 0])
    alpha = circ_ang2rad(PDs);
    alpha = alpha(~isnan(alpha));
    [pval, ~] = circ_rtest(alpha);
    if pval < 0.001
        fprintf('Within anisotropic zone: N=%i, p=%.3e (Rayleigh test) \n', numel(PDs), pval);
    else
        fprintf('Within anisotropic zone: N=%i, p=%.3f (Rayleigh test) \n', numel(PDs), pval);
    end
    angDev_CH(1,mouseId) = concent_param(alpha, 0);
    
    ind = AnisotropyInd>1 & subMask;
    inFlag_fS1 = inpolygon(CenterXY_CH(ind,1), CenterXY_CH(ind,2), Suc.xySensory(:,1), Suc.xySensory(:,2));
    inFlag_CFA = inpolygon(CenterXY_CH(ind,1), CenterXY_CH(ind,2), Suc.xyMotor(:,1), Suc.xyMotor(:,2));
    Nin = sum(inFlag_fS1 | inFlag_CFA);
    fprintf('Within fS1/CFA: N=%i/%i \n', Nin, numel(PDs))
    Fract_in_fS1_CFA_CH(1,mouseId) = Nin/numel(PDs);
    
    PDs = PDs_CH(AnisotropyInd>1 & subMask & CenterXY_CH(:,1)<-1.75);
    subplot(4,Nmice,mouseId+2*Nmice)
    plotPDs(PDs, gca, 1, 0, [0 0.6 1])
    fprintf('Lateral neurons: N=%i \n', numel(PDs))
    alpha1 = circ_ang2rad(PDs);
    alpha1 = alpha1(~isnan(alpha1));

    PDs = PDs_CH(AnisotropyInd>1 & subMask & CenterXY_CH(:,1)>=-1.75);
    subplot(4,Nmice,mouseId+3*Nmice)
    plotPDs(PDs, gca, 1, 0, [0.97 0 1])
    fprintf('Medial neurons: N=%i \n', numel(PDs))
    alpha2 = circ_ang2rad(PDs);
    alpha2 = alpha2(~isnan(alpha2));

    [pval, ~] = circ_wwtest(alpha1, alpha2);
    fprintf('Wason-Williams test for equal means (Lateral vs. Medial): p=%.2e', pval);
    
    angMean_CH(:,mouseId) = [circ_mean(alpha1);circ_mean(alpha2)];

    fprintf('\n\n\n');


end


%% Figure S1C


angDev_EE=mod(circ_rad2ang(angDev_EE), 360);
angDev_CH=mod(circ_rad2ang(angDev_CH), 360);


figure('Color', 'w', 'Position',[1400 600 250 300], 'Name', 'Supp Figure 1C')
plot((rand(1,5)-0.5)+2, angDev_EE, 'o', 'MarkerSize',8, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor','none')
hold on, plot((rand(1,4)-0.5)+4, angDev_CH, 'o', 'MarkerSize',8, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor','none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [1 5], 'YLim', [50 90], 'YTick', 0:0.25:1, 'XTick', [2 4], 'XTickLabel', {'EE mice'; 'CH mice'})
ylabel('Average angular deviation from the mean (deg)')
p = permutationTest(angDev_EE, angDev_CH);
pTxt = sprintf('p=%.3f', p);
text(3,85,pTxt);


%% Figure S1D


figure('Color', 'w', 'Position',[1400 400 250 300], 'Name', 'Supp Figure 1D')
plot((rand(1,5)-0.5)+2, Fract_in_fS1_CFA_EE, 'o', 'MarkerSize',8, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor','none')
hold on, plot((rand(1,4)-0.5)+4, Fract_in_fS1_CFA_CH, 'o', 'MarkerSize',8, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor','none')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [1 5], 'YLim', [0.5 1], 'YTick', 0:0.25:1, 'XTick', [2 4], 'XTickLabel', {'EE mice'; 'CH mice'})
ylabel('Fraction of cells within fS1/CFA boundaries')
p = permutationTest(Fract_in_fS1_CFA_EE, Fract_in_fS1_CFA_CH);
pTxt = sprintf('p=%.3f', p);
text(3,0.95,pTxt);


%% Figure S1E

angMean_EE =  mod(circ_rad2ang(angMean_EE), 360);
angMean_CH = mod(circ_rad2ang(angMean_CH), 360);

figure('Color', 'w', 'Position',[1400 200 250 300], 'Name', 'Supp Figure 1E')
for k=1:size(angMean_EE,2)
    plot([1 2], [angMean_EE(1,k) angMean_EE(2,k)], 'Marker','none', 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
    hold on
end
for k=1:size(angMean_CH,2)
    plot([3 4], [angMean_CH(1,k) angMean_CH(2,k)], 'Marker','none', 'Color', [0.7 0.7 0.7], 'LineWidth',1)
end
plot(ones(1,5), angMean_EE(1,:), 'o', 'MarkerEdgeColor',[0 0.6 1], 'MarkerFaceColor',[0 0.6 1], 'LineWidth',1, 'MarkerSize', 6)
plot(2*ones(1,5), angMean_EE(2,:), 'o', 'MarkerEdgeColor',[0.97 0 1], 'MarkerFaceColor',[0.97 0 1], 'LineWidth',1, 'MarkerSize', 6)
plot(3*ones(1,4), angMean_CH(1,:), 'o', 'MarkerEdgeColor',[0 0.6 1], 'MarkerFaceColor','w', 'LineWidth',1, 'MarkerSize', 6)
plot(4*ones(1,4), angMean_CH(2,:), 'o', 'MarkerEdgeColor',[0.97 0 1], 'MarkerFaceColor','w', 'LineWidth',1, 'MarkerSize', 6)
hold off
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Xlim', [0 5], 'XTick', [1.5 3.5], 'XTickLabel', {'EE mice'; 'CH mice'}, 'YLim', [150 400])
ylabel('Mean preferred direction (deg)')

[~,p] = ttest(angMean_EE(1,:)-angMean_EE(2,:));
fprintf('Mean preferred directions, EE mice (Lateral vs. Medial): p=%.2f (t-test) \n', p);
[~,p] = ttest(angMean_CH(1,:)-angMean_CH(2,:));
fprintf('Mean preferred directions, CH mice (Lateral vs. Medial): p=%.2e (t-test) \n', p);
p = permutationTest(angMean_EE(1,:), angMean_CH(1,:));
fprintf('Mean preferred directions, Lateral neurons (EE vs. CH): p=%.3f (permutation test) \n', p);
p = permutationTest(angMean_EE(2,:), angMean_CH(2,:));
fprintf('Mean preferred directions, Medial neurons (EE vs. CH): p=%.3f (permutation test) \n', p);