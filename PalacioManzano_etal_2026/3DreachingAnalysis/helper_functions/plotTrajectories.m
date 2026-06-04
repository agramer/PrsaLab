function plotTrajectories(xyzAll, hAx, col, normFlag, trialInds)


hold(hAx, 'On');

Ntrials = size(xyzAll{1,1},1);

if nargin<5 || isempty(trialInds)
    trialInds = 1:Ntrials;
end

for tr=trialInds
    plot3(hAx, xyzAll{1,1}(tr,:), xyzAll{1,2}(tr,:), xyzAll{1,3}(tr,:), 'Color', col, 'LineWidth', 1)
end
axis equal
if normFlag
    set(hAx, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'YDir', 'Reverse', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5, ...
        'XLim', [-10 25], 'YLim', [-13 13], 'ZLim', [-5 20], 'View', [-130 25], 'XTick', -10:5:25, 'YTick', -10:5:10, 'ZTick', -5:5:20)
    set(hAx.XAxis, 'Visible', 'off');
    set(hAx.YAxis, 'Visible', 'off');
    set(hAx.ZAxis, 'Visible', 'off');
    
    hold(hAx, 'On')
    line(hAx, [-10 0], [13 13], [-5 -5], 'Color', [0 0.8 0], 'LineWidth', 2)
    line(hAx, [-10 -10], [13 3], [-5 -5], 'Color', [0.8 0 0], 'LineWidth', 2)
    line(hAx, [-10 -10], [13 13], [-5 5], 'Color', [0 0 0.8], 'LineWidth', 2)
    text(hAx, -22,4,7, '1 cm', 'FontSize',12,'FontWeight','bold')
else
    set(hAx, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'YDir', 'Reverse', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5,...
        'XLim', [-25 5], 'YLim', [-11 9], 'ZLim', [-25 10], 'View', [-50 30])
end
hold(hAx, 'Off');


