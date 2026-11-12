function [mu_deg, kappa, r_obs, p_value] = fitVonMises(angles_deg, values, numPermutations, plotFlag)
% Fit a von Mises distribution to angular data in degrees and plot results in both Cartesian and polar coordinates.
%
% INPUTS:
%   angles_deg      - vector of angles in degrees
%   values          - corresponding values (same size)
%   numPermutations - (optional) number of permutations for significance test (default = 1000)
%
% OUTPUTS:
%   mu_deg  - estimated mean direction (degrees)
%   kappa   - concentration parameter
%   r_obs   - Pearson correlation between fit and data
%   p_value - significance of the fit (permutation test)

% See Fisher, N. I. (1995). Statistical Analysis of Circular Data. Cambridge University Press.
% Computations of Von Mises distribution parameters and pdf are available
% in the CirStat Matlab toolbox by Berens

    if nargin < 3
        numPermutations = 1000;
    end

    % Ensure column vectors
    angles_deg = angles_deg(:);
    values = values(:);

    % Convert degrees to radians
    angles_rad = deg2rad(angles_deg);

    % Normalize values
    values = values / sum(values);

    % Compute weighted vector
    C = sum(values .* cos(angles_rad));
    S = sum(values .* sin(angles_rad));
    mu_rad = atan2(S, C);
    if mu_rad < 0, mu_rad = mu_rad + 2*pi; end
    mu_deg = rad2deg(mu_rad);

    R = sqrt(C^2 + S^2);

    % Estimate kappa
    if R < 0.53
        kappa = 2 * R + R^3 + (5 * R^5) / 6;
    elseif R < 0.85
        kappa = -0.4 + 1.39 * R + 0.43 / (1 - R);
    else
        kappa = 1 / (R^3 - 4 * R^2 + 3 * R);
    end
    if kappa < 1e-6, kappa = 0; end

    % Von Mises PDF
    fitted_pdf = exp(kappa * cos(angles_rad - mu_rad)) / (2*pi * besseli(0, kappa));
    fitted_pdf = fitted_pdf * max(values) / max(fitted_pdf);  % scale to match data peak

    % Observed correlation
    r_obs = corr(values, fitted_pdf, 'Type', 'Pearson');

    % Permutation test
    r_null = zeros(numPermutations, 1);
    for i = 1:numPermutations
        shuffled = values(randperm(length(values)));
        r_null(i) = corr(shuffled, fitted_pdf, 'Type', 'Pearson');
    end
    p_value = mean(r_null >= r_obs);
    
    if plotFlag
        % Generate von Mises PDF
        theta = linspace(0, 2*pi-pi/4, 200)';
        
        pdf = exp(kappa * cos(theta - mu_rad)) / (2*pi * besseli(0, kappa));
        % Scale PDF to match the peak of input data for visual comparison
        scaled_pdf = pdf * max(values) / max(pdf);

        figure('Color', 'w')
        % Plot: Cartesian
        subplot(1,2,1);
        plot(angles_deg, values, '-ko', 'DisplayName', 'Data'); hold on;
        plot(rad2deg(theta), scaled_pdf, 'r-', 'LineWidth', 2, 'DisplayName', 'Von Mises Fit');
        xlabel('Angle (°)');
        ylabel('Normalized value');
        title(sprintf('Cartesian Fit: \\mu = %.1f°, \\kappa = %.2f, r = %.2f, p = %.3f', ...
                      mu_deg, kappa, r_obs, p_value));
        legend('Location', 'best');
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03]);
    

        % Generate von Mises PDF
        theta = linspace(0, 2*pi, 200)';
        
        pdf = exp(kappa * cos(theta - mu_rad)) / (2*pi * besseli(0, kappa));
        % Scale PDF to match the peak of input data for visual comparison
        scaled_pdf = pdf * max(values) / max(pdf);

        % Plot: Polar
        subplot(1,2,2);
        polarplot([angles_rad; angles_rad(1)], [values;values(1)], '-ko'); hold on;
        polarplot([theta; theta(1)], [scaled_pdf; scaled_pdf(1)], 'r-', 'LineWidth', 2);
        title('Polar Plot of Von Mises Fit');
    end

end