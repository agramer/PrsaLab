function plotDirTune(trigData, cellNum, hFig, col)


plotMap = [6 3 2 1 4 7 8 9];
angles = [0:pi/4:2*pi-pi/4]';

UpLim = [];
LoLim = [];

figure(hFig);
for dir=1:8
    subplot(3,3,plotMap(dir))
    hold on

    idx = (trigData.Info(:,4)==dir & (trigData.Info(:,3)==5 | trigData.Info(:,3)==4 | trigData.Info(:,3)==3)...
        & trigData.Info(:,2)==1);   % outward mvt
    if sum(idx)>0
        plot(trigData.periTime,mean(trigData.Traces{cellNum,1}(idx,:)), 'Color', col, 'LineWidth', 2)

        % Create closed polygon coordinates
        x_shade = [trigData.periTime, fliplr(trigData.periTime)];
        y_shade = [mean(trigData.Traces{cellNum,1}(idx,:))-std(trigData.Traces{cellNum,1}(idx,:)), fliplr(mean(trigData.Traces{cellNum,1}(idx,:))+std(trigData.Traces{cellNum,1}(idx,:)))];
        fill(x_shade, y_shade, col, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
    hold off
    set(gca, 'Box', 'off', 'XLim', [-1 1.5], 'XColor', 'none', 'YColor', 'none')
    
    switch dir
        case 5
            ylabel('df/f')
        case 7
            xlabel('Time (s)')
    end
    if sum(idx)>0
        ylims = get(gca, 'YLim');
        UpLim = [UpLim ylims(2)];
        LoLim = [LoLim ylims(1)];
    end
    
end

[X,Y] = pol2cart(angles,trigData.Stats{cellNum,1}.PeakDff(:,1));
X = [X;X(1)];
Y = [Y;Y(1)];

subplot(3,3,5)
plot(X,Y,'o-k', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 4, 'LineWidth', 1)
line([-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], [0 0], 'LineStyle', ':', 'LineWidth', 0.5, 'Color', 'k')
line([0 0], [-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], 'LineStyle', ':', 'LineWidth', 0.5, 'Color', 'k')
line([-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], [-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], 'LineStyle', ':', 'LineWidth', 0.5, 'Color', 'k')
line([-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], [max(abs([X;Y]))+0.05 -max(abs([X;Y]))-0.05], 'LineStyle', ':', 'LineWidth', 0.5, 'Color', 'k')

axis equal 
set(gca, 'Box', 'off', 'XLim', [-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], 'YLim', [-max(abs([X;Y]))-0.05 max(abs([X;Y]))+0.05], 'XColor', 'none', 'YColor', 'none')


UpLim=max(UpLim);
LoLim=min(LoLim);

for dir=1:8
    figure(hFig)
    subplot(3,3,plotMap(dir))
    line([0 0], [LoLim-0.1 UpLim+0.1], 'LineStyle', ':', 'LineWidth', 1, 'Color', 'k')
    set(gca, 'YLim', [LoLim-0.1 UpLim+0.1])
end





    