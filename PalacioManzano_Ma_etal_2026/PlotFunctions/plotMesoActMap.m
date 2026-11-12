function plotMesoActMap(ANAframe, Suc, hAx)

plDIM = [4 4];
imDIM = [201 201];
pRES = plDIM./(imDIM-1);
Xref = -plDIM(1):pRES(1):0;
Yref = -plDIM(2)/2:pRES(2):plDIM(2)/2;

hold(hAx, 'on')

Cmap = fire(256);

winStep = pRES(1);
plot3(hAx, Suc.xySensory(:,1)+winStep/2, Suc.xySensory(:,2)+winStep/2, 3*ones(1,size(Suc.xySensory,1)), 'w--')
plot3(hAx, Suc.xyMotor(:,1)+winStep/2, Suc.xyMotor(:,2)+winStep/2, 3*ones(1,size(Suc.xyMotor,1)), 'w--')
surf(hAx, Xref, Yref, ANAframe);
shading flat
hold off
set(hAx, 'XLim', [-3.5 Xref(end)]+winStep/2, 'YLim', [Yref(1) Yref(end)]+winStep/2, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XTick', (Xref(1):0.5:Xref(end))+winStep/2, ...
    'YTick', (Yref(1):0.5:Yref(end))+winStep/2, 'XTickLabel',  Xref(1):0.5:Xref(end),  'YTickLabel', Yref(1):0.5:Yref(end), 'XGrid', 'on', 'YGrid', 'on', 'Box', 'off', 'CLim', [0 max(ANAframe(:))])
xlabel('Medio-lateral position (mm)')
ylabel('Antero-posterior position (mm)')
view(0,90)
axis square
colormap(gca, Cmap);

hold(hAx, 'off')