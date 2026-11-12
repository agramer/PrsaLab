function stats = plotMechEfficiency(dat1, dat2, doPlot, removeOutliers)
% plotMechEfficiency
% Fits linear regression dat2 ~ dat1.
% Ignores NaNs pairwise.
% Optionally removes outliers before fitting.
%
% INPUT:
%   dat1, dat2       : Nsamples x 1 vector
%   doPlot           : (optional) true/false, default = true
%   removeOutliers   : (optional) true/false, default = false
%
% OUTPUT:
%   stats : structure with regression results 

    if nargin < 3 || isempty(doPlot)
        doPlot = true;
    end
    if nargin < 4 || isempty(removeOutliers)
        removeOutliers = false;
    end

    if ~isequal(size(dat1), size(dat2))
        error('dat1 and dat2 must have the same size.');
    end


    stats = struct();
    stats.session = struct('slope', NaN, 'intercept', NaN, ...
        'R2', NaN, 'p', NaN, 'N', 0, 'Nout', 0);

    if doPlot
        figure('Name','dat1 vs dat2','Color','w');
    end

    
    x = dat1;
    y = dat2;

    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    outMask = false(size(x));

    if numel(x) >= 3 && removeOutliers
        % Remove outliers based on residuals from an initial fit
        mdl0 = fitlm(x, y);
        outMask = mdl0.Diagnostics.CooksDistance > (4/numel(x));
        xFit = x(~outMask);
        yFit = y(~outMask);
    else
        xFit = x;
        yFit = y;
    end

    if numel(xFit) >= 2
        mdl = fitlm(xFit, yFit);

        slope = mdl.Coefficients.Estimate(2);
        intercept = mdl.Coefficients.Estimate(1);
        R2 = mdl.Rsquared.Ordinary;
        pval = mdl.Coefficients.pValue(2);
    else
        slope = NaN;
        intercept = NaN;
        R2 = NaN;
        pval = NaN;
    end

    % Store stats
    stats.slope = slope;
    stats.intercept = intercept;
    stats.R2 = R2;
    stats.p = pval;
    stats.N = numel(xFit);
    stats.Nout = sum(outMask);
    stats.MeanX = mean(xFit);
    stats.xData = xFit;
    stats.yData = yFit;
    

    % Plot if requested
    if doPlot
        % Plot included points
        scatter(xFit, yFit, 25, 'filled');
        hold on

        % Plot excluded outliers (if any)
        if removeOutliers && any(outMask)
            scatter(x(outMask), y(outMask), 40, 'o', 'LineWidth', 1.5);
        end

        if numel(xFit) >= 2
            xx = linspace(min(xFit), max(xFit), 100);
            yy = intercept + slope * xx;
            plot(xx, yy, 'LineWidth', 2);

            title(sprintf('slope=%.3g, R^2=%.3g, p=%.3g, N=%d (out=%d)', ...
                slope, R2, pval, numel(xFit), sum(outMask)));
        else
            title(sprintf('not enough data'));
        end

        xlabel('dat1');
        ylabel('dat2');
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XGrid', 'off', 'YGrid', 'off')
    end
    
end