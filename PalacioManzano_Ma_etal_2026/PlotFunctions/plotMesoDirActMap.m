function plotMesoDirActMap(ANAframes_DIR, Suc, mouseID, hFig)

Ndirs = 8;

plDIM = [4 4];
imDIM = [201 201];
pRES = plDIM./(imDIM-1);
Xref = -plDIM(1):pRES(1):0;
Yref = -plDIM(2)/2:pRES(2):plDIM(2)/2;
winStep = pRES(1);

Cmap = fire(256);
figure(hFig);

indSP = [6 3 2 1 4 7 8 9];

for dr=1:Ndirs
    subplot(3,3,indSP(dr))
    hold on

    plot3(Suc.xySensory(:,1)+winStep/2, Suc.xySensory(:,2)+winStep/2, 3*ones(1,size(Suc.xySensory,1)), 'w--')
    plot3(Suc.xyMotor(:,1)+winStep/2, Suc.xyMotor(:,2)+winStep/2, 3*ones(1,size(Suc.xyMotor,1)), 'w--')
    surf(Xref, Yref, ANAframes_DIR{dr,1}(:,:,mouseID));
    shading flat
    hold off
    set(gca, 'XLim', [-3.5 Xref(end)]+winStep/2, 'YLim', [Yref(1) Yref(end)]+winStep/2, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XTick', (Xref(1):0.5:Xref(end))+winStep/2, ...
        'YTick', (Yref(1):0.5:Yref(end))+winStep/2, 'XTickLabel',  Xref(1):0.5:Xref(end),  'YTickLabel', Yref(1):0.5:Yref(end), 'XGrid', 'on', 'YGrid', 'on', 'Box', 'off', 'CLim', [0 1])
    xlabel('Medio-lateral position (mm)')
    ylabel('Antero-posterior position (mm)')
    view(0,90)
    axis square

    colormap(Cmap)
end