function AnisotropyInd = plotTopography(PDs, CenterXY, Suc, meanFlag, hAx)

MLlims = [-4 0.5];
APlims = [-2 2];

winLen = 0.8;
winStep = 0.05;
cellNumThr = 10;

MLcoords = MLlims(1):winStep:MLlims(2);
APcoords = APlims(1):winStep:APlims(2);
pRayleigh = zeros(numel(MLcoords), numel(APcoords));
cMean = -ones(numel(MLcoords), numel(APcoords));


for ml=1:numel(MLcoords)
    for ap=1:numel(APcoords)
        dist = sqrt( (CenterXY(:,1)-MLcoords(ml)).^2 + (CenterXY(:,2)-APcoords(ap)).^2 );
        
        idxCells = dist<=winLen/2; 
        alpha = circ_ang2rad(PDs(idxCells,1));

        if sum(~isnan(alpha))>=cellNumThr
            [pVal, ~] = circ_rtest(alpha(~isnan(alpha)));
            if pVal<0.01
                pRayleigh(ml,ap) = 3;
                mVal = circ_mean(alpha(~isnan(alpha)));
                mVal = mod(circ_rad2ang(mVal), 360);
            elseif pVal<0.05
                pRayleigh(ml,ap) = 2;
                mVal = circ_mean(alpha(~isnan(alpha)));
                mVal = mod(circ_rad2ang(mVal), 360);
            else
                pRayleigh(ml,ap) = 1;
                mVal = -1;
            end
            cMean(ml,ap) = mVal;
        end
        

    end
end

Ncells = size(PDs,1);
AnisotropyInd = zeros(Ncells, 1);
for k=1:Ncells
    [~,indML] = min(abs(CenterXY(k,1)-MLcoords));
    [~,indAP] = min(abs(CenterXY(k,2)-APcoords));
    AnisotropyInd(k, 1) = pRayleigh(indML,indAP);
end

if meanFlag
    data = cMean;
    cmap=hsv(360);
    cmap=[ones(1,3);cmap];
    c_Lim = [0 360];
else
    data = pRayleigh;
    cmap=parula(3);
    cmap=[ones(1,3);cmap];
    c_Lim = [0 3];
end

data = [data data(:,end)];
data = [data; data(end,:)];

surf(hAx, [MLcoords MLcoords(end)], [APcoords APcoords(end)], data')
colormap(cmap)
shading flat

hold(hAx, 'on')
plot3(hAx, Suc.xySensory(:,1)+winStep/2, Suc.xySensory(:,2)+winStep/2, 3*ones(1,size(Suc.xySensory,1)), 'k--')
plot3(hAx, Suc.xyMotor(:,1)+winStep/2, Suc.xyMotor(:,2)+winStep/2, 3*ones(1,size(Suc.xyMotor,1)), 'k--')
hold(hAx, 'off')
view(0,90)
set(hAx, 'XLim', [-3.5 0]+winStep/2, 'YLim', APlims+winStep/2, 'XTick', (MLlims(1):0.5:MLlims(2))+winStep/2, 'XTickLabel', (MLlims(1):0.5:MLlims(2)), ...
    'YTick', (APlims(1):0.5:APlims(2))+winStep/2, 'YTickLabel', (APlims(1):0.5:APlims(2)), 'TickDir', 'out', 'TickLength', [0.03 0.03], ...
    'Clim', c_Lim, 'XGrid', 'on', 'YGrid', 'on')
xlabel('Medio-lateral position (mm)')
ylabel('Antero-posterior position (mm)')

