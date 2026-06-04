function regStats = plotTrialKinematics(hAx, varVals, col, varNames, varIDa, varIDb)


hold(hAx, 'On');

Ntrials = size(varVals,1);

plot(varVals(:, varIDa), varVals(:, varIDb), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none')
[B,~,~,~,STATS] = regress(varVals(:, varIDb),[ones(Ntrials,1) varVals(:,varIDa)]);
line([min(varVals(:,varIDa)) max(varVals(:,varIDa))], B(2)*[min(varVals(:,varIDa)) max(varVals(:,varIDa))]+B(1), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2)
hold off
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03])
ylabel(varNames{varIDb})
xlabel(varNames{varIDa})
title(['r2=' num2str(STATS(1)) ' p=' num2str(STATS(3))])
axis square
hold(hAx, 'Off');

regStats = [min(varVals(:,varIDa)) max(varVals(:,varIDa)) B' STATS([1 3])];