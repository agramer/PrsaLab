function plotMapAvgAct(maps1, maps2)

plDIM = [4 4];
imDIM = [201 201];
pRES = plDIM./(imDIM-1);
Xref = -plDIM(1):pRES(1):0;
Yref = -plDIM(2)/2:pRES(2):plDIM(2)/2;

figure('Color', 'w', 'Position', [210 300 800 600], 'Name', 'Figure 2E');
datML_ctr(:,:) = mean(maps1,1);
datAP_ctr(:,:) = mean(maps1,2);
datML_enr(:,:) = mean(maps2,1);
datAP_enr(:,:) = mean(maps2,2);

subplot(4,2,5:2:7)
hold on
plot(Xref, datML_ctr, 'Color', [0.7 0.7 0.7])
plot(Xref, datML_enr, 'Color', [1 0.7 0.7])
plot(Xref, mean(datML_ctr,2), 'Color', 'k', 'LineWidth', 2)
plot(Xref, mean(datML_enr,2), 'Color', 'r', 'LineWidth', 2)
hold off
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 Xref(end)], 'YLim', [0 0.52])
xlabel('Medio lateral position (mm)')
ylabel('Average activity')

subplot(4,2,6:2:8)
hold on
plot(Yref, datAP_ctr, 'Color', [0.7 0.7 0.7])
plot(Yref, datAP_enr, 'Color', [1 0.7 0.7])
plot(Yref, mean(datAP_ctr,2), 'Color', 'k', 'LineWidth', 2)
plot(Yref, mean(datAP_enr,2), 'Color', 'r', 'LineWidth', 2)
hold off
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [Yref(1) Yref(end)], 'YLim', [0 0.52])
xlabel('Antero posterior position (mm)')
ylabel('Average activity')


NsessCtr = size(maps1,3);
NsessEnr = size(maps2,3);

maxXYctr = zeros(NsessCtr,2);
maxXYenr = zeros(NsessEnr,2);
areaCtr = zeros(NsessCtr,1);
areaEnr = zeros(NsessEnr,1);

for k=1:NsessCtr
    clear dat
    dat(:,:) = maps1(:,:,k); 
    [~,ind] = max(dat,[],'all');
    [yi,xi] = ind2sub(size(dat),ind);       % xi is column (2nd index), yi is row (1st index)
    maxXYctr(k,:) = [Xref(xi) Yref(yi)];    
    areaCtr(k,1) = plDIM(1)*plDIM(2)*sum(dat(:)>0)/(size(dat,1)*size(dat,2));
end

for k=1:NsessEnr
    clear dat
    dat(:,:) = maps2(:,:,k); 
    [~,ind] = max(dat,[],'all');
    [yi,xi] = ind2sub(size(dat),ind);
    maxXYenr(k,:) = [Xref(xi) Yref(yi)];    % xi is column (2nd index), yi is row (1st index)
    areaEnr(k,1) = plDIM(1)*plDIM(2)*sum(dat(:)>0)/(size(dat,1)*size(dat,2));
end

subplot(4,2,1)
hold on
plot(maxXYctr(:,1), rand(1,NsessCtr)-0.5+1.5, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.7 0.7 0.7])
plot(maxXYenr(:,1), rand(1,NsessEnr)-0.5, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0.7 0.7])
plot(mean(maxXYctr(:,1)), 1.5, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5)
plot(mean(maxXYenr(:,1)), 0, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none', 'LineWidth', 1.5)
hold off
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'YTick', [0 2], 'YTickLabel', {'EE'; 'CH'},  'XLim', [-3.5 Xref(end)], 'YLim', [-1 2.5], 'Box', 'off')
xlabel('Medio lateral position (mm)')

subplot(4,2,2)
hold on
plot(maxXYctr(:,2), rand(1,NsessCtr)-0.5+1.5, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.7 0.7 0.7])
plot(maxXYenr(:,2), rand(1,NsessEnr)-0.5, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0.7 0.7])
plot(mean(maxXYctr(:,2)), 1.5, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5)
plot(mean(maxXYenr(:,2)), 0, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none', 'LineWidth', 1.5)
hold off
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'YTick', [0 2], 'YTickLabel', {'EE'; 'CH'},  'XLim', [Yref(1) Yref(end)], 'YLim', [-1 2.5], 'Box', 'off')
xlabel('Antero posterior position (mm)')

[~,p] = ttest2(maxXYctr(:,1), maxXYenr(:,1));
disp(['Max ML location CH vs. EE (two-sample t test): p=' num2str(p)])

[~,p] = ttest2(maxXYctr(:,2), maxXYenr(:,2));
disp(['Max AP location CH vs. EE (two-sample t test): p=' num2str(p)])

[~,p] = ttest2(areaCtr, areaEnr);
disp(['Activation area CH vs. EE (two-sample t test): p=' num2str(p)])


subplot(4,2,3)
AUC = runROC(Xref, datML_ctr', datML_enr', 0.05, 999, gca);
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-3.5 Xref(end)], 'YLim', [0 1], 'Box', 'off')
xlabel('Medio lateral position (mm)')
subplot(4,2,4)
AUC = runROC(Yref, datAP_ctr', datAP_enr', 0.05, 999, gca);
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [Yref(1) Yref(end)], 'YLim', [0 1], 'Box', 'off')
xlabel('Medio lateral position (mm)')