function plotMetabData(data, idxEE, idxCH, yLims, yLbl, hAx)

nEE = sum(idxEE);
nCH = sum(idxCH);

hold(hAx, 'on');
plot(hAx, rand(1,nCH)-0.5, data(idxCH), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none')
plot(hAx, rand(1,nEE)-0.5+2, data(idxEE), 'o', 'MarkerFaceColor', [1 0.7 0.7], 'MarkerEdgeColor', 'none')
errorbar(hAx, 0.75, mean(data(idxCH), 'omitnan'), std(data(idxCH), 'omitnan'), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [0 0 0])
errorbar(hAx, 2.75, mean(data(idxEE), 'omitnan'), std(data(idxEE), 'omitnan'), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [1 0 0])
set(hAx, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-1 3], 'YLim', yLims, 'XTick', [0 2], 'XTickLabel', {'CH'; 'EE'})
ylabel(yLbl)
hold(hAx, 'off');



