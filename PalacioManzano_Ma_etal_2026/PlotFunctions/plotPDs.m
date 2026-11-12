function plotPDs(PDs, hAx, histFlag, meanFlag, col, R)
% R (optional): output of circ_vmcomponents for these PDs; the von Mises
% mixture of the selected model (R.components) is drawn on top of the
% histogram, scaled to the same expected counts per bin.


PDs=mod(PDs, 360);
[x,y]=pol2cart(deg2rad(PDs),1);

hold(hAx, 'on')
rgbA = applyAlpha(col, 0.2);
hP = plot(hAx, [zeros(size(x)) x]',[zeros(size(x)) y]', 'Color', rgbA);

if meanFlag
    PDmean = mod(circ_rad2ang(circ_mean(circ_ang2rad(PDs))), 360);
    [x,y]=pol2cart(deg2rad(PDmean),1);
    plot(hAx, [0 x], [0 y], 'Color', 'r', 'LineWidth', 2, 'LineStyle', ':')
end
axis square
axis equal
set(hAx, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XLim', [-1.2 1.2], 'YLim', [-1.2 1.2])
hold(hAx, 'off')


if histFlag
    hold(hAx, 'on')
    Nbins = 12;
    histRes = 360/Nbins;
    hVals = hist(PDs+360, 360:histRes:360+360-histRes);
    
    circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];         % Circle Function For Angles In Radians
    circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];       % Circle Function For Angles In Degrees
    N = 360*2;                                                                      % Number Of Points In Complete Circle
    
    piFrac = deg2rad(histRes/2);
    for bin=1:Nbins
        r_angl = linspace((2*bin-3)*piFrac, (2*bin-1)*piFrac, N);                       % Angle Defining Arc Segment (radians)
        radius = hVals(bin)/max(hVals);                                                 % Arc Radius
        xy_r = circr(radius,r_angl);                                                    % Matrix (2xN) Of (x,y) Coordinates
        plot(hAx, xy_r(1,:), xy_r(2,:), 'Color', col, 'LineWidth', 1.5)                                   % Draw an arc
        [x,y] = pol2cart((2*bin-3)*piFrac,radius);
        plot(hAx, [0 x], [0 y], 'Color', col, 'LineWidth', 1.5)
        [x,y] = pol2cart((2*bin-1)*piFrac,radius);
        plot(hAx, [0 x], [0 y], 'Color', col, 'LineWidth', 1.5)
    end
    hold(hAx, 'off')

end


if nargin > 5 && ~isempty(R)
    hold(hAx, 'on')
    Nbins = 12;
    histRes = 360/Nbins;
    hVals = hist(PDs+360, 360:histRes:360+360-histRes);                             % same bins as the histogram (scaling)
    C = R.components;
    loc = C.location;
    if all(loc <= 2*pi) && R.meanDir <= 2*pi, loc = rad2deg(loc); end               % R computed in radians
    vm = @(t, j) C.weight(j) * exp(C.kappa(j)*(cos(t - deg2rad(loc(j))) - 1)) / (2*pi*besseli(0, C.kappa(j), 1));
    t = linspace(0, 2*pi, 721);
    f = zeros(size(t));
    for j = 1:height(C), f = f + vm(t, j); end                                     % mixture of the selected components
    sc = deg2rad(histRes) * sum(~isnan(PDs)) / max(hVals);                          % density -> expected count per bin, normalized as the histogram
    [x,y] = pol2cart(t, f*sc);
    plot(hAx, x, y, 'k:', 'LineWidth', 1.5)
    for j = 1:height(C)                                                             % markers at the component means
        fj = 0;  for k = 1:height(C), fj = fj + vm(deg2rad(loc(j)), k); end
        [x,y] = pol2cart(deg2rad(loc(j)), fj*sc);
        plot(hAx, x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
    end
    hold(hAx, 'off')

    % text summary of the mixture fit
    fprintf('   von Mises mixture (n = %d neurons): %d component(s) selected by BIC\n', sum(~isnan(PDs)), height(C));
    for j = 1:height(C)
        fprintf('   component %d: %.0f%% of neurons, mean %.0f deg, kappa %.2f\n', ...
            j, 100*C.weight(j), mod(loc(j), 360), C.kappa(j));
    end
    dB = R.BIC - R.BIC(height(C));
    fprintf('   BIC:');
    for K = 1:numel(R.BIC), fprintf('  K = %d: %.1f (dBIC %+.1f)', K, R.BIC(K), dB(K)); end
    fprintf('\n   (dBIC > 0: model worse than the selected one)\n');
end