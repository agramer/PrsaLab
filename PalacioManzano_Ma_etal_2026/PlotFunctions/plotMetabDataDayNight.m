function plotMetabDataDayNight(dataDay, dataNight, idxEEday, idxCHday, idxEEnight, idxCHnight, yLims, yLbl, hAx)

nEEday = sum(idxEEday);
nCHday = sum(idxCHday);

nEEnight = sum(idxEEnight);
nCHnight = sum(idxCHnight);

hold(hAx, 'on');
rectangle('Position', [-1 yLims(1) 4.5 yLims(2)-yLims(1)], 'FaceColor', [0.98 0.96 0.82], 'EdgeColor', 'none');
rectangle('Position', [3.5 yLims(1) 4.5 yLims(2)-yLims(1)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');


plot(hAx, rand(1,nCHday)-0.5, dataDay(idxCHday), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none')
plot(hAx, rand(1,nEEday)-0.5+2, dataDay(idxEEday), 'o', 'MarkerFaceColor', [1 0.7 0.7], 'MarkerEdgeColor', 'none')
errorbar(hAx, 0.75, mean(dataDay(idxCHday)), std(dataDay(idxCHday)), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [0 0 0])
errorbar(hAx, 2.75, mean(dataDay(idxEEday)), std(dataDay(idxEEday)), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [1 0 0])

plot(hAx, rand(1,nCHnight)-0.5+5, dataNight(idxCHnight), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none')
plot(hAx, rand(1,nEEnight)-0.5+7, dataNight(idxEEnight), 'o', 'MarkerFaceColor', [1 0.7 0.7], 'MarkerEdgeColor', 'none')
errorbar(hAx, 5.75, mean(dataNight(idxCHnight)), std(dataNight(idxCHnight)), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [0 0 0])
errorbar(hAx, 7.75, mean(dataNight(idxEEnight)), std(dataNight(idxEEnight)), 'CapSize', 0,  'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [1 0 0])

set(hAx, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-1 8], 'YLim', yLims, 'XTick', [0 2 5 7], 'XTickLabel', {'CH'; 'EE'; 'CH'; 'EE'})
ylabel(yLbl)
hold(hAx, 'off');


