function [p_vals, obs_stats] = cluster_permutation_test(maps1, maps2, n_perm, alpha, Suc)
% CLUSTER_PERMUTATION_TEST_WITH_PLOTS 
% Performs a non-parametric permutation test for small sample sizes (N=9 vs N=7).
%
% DEFINITION OF THE T-VALUE:
% The t-value (t-statistic) represents the "Signal-to-Noise" ratio for each pixel.
%   t = (Mean_Group1 - Mean_Group2) / Standard_Error
%
% THE THREE TEST STATISTICS:
% 1. MAX PIXEL T-STAT: Highest absolute t-value at any single pixel.
% 2. MAX CLUSTER SIZE: Total pixels in the largest connected group exceeding t_thresh.
% 3. MAX CLUSTER MASS: Sum of absolute t-values within the largest cluster.
%
% INPUTS:
%   maps1:    [201x201x9] Double - Control group maps
%   maps2:    [201x201x7] Double - Enriched group maps
%   n_perm:   (Optional) Number of shuffles. Default: 5000.
%   t_thresh: (Optional) T-cutoff. Default: 2.14 (p<0.05 for df=14, df = N1+N2-2).

    % --- Handle Default Values for N=9 vs N=7 ---
    if nargin < 3 || isempty(n_perm)
        n_perm = 5000;
    end
    if nargin < 4 || isempty(alpha)
        alpha = 0.01; % Significance level
    end


    [nx, ny, n1] = size(maps1);
    [~, ~, n2] = size(maps2);
    all_maps = cat(3, maps1, maps2);
    total_n = n1 + n2;

    % calculate the t-threshold for the desired significance of the
    % t-statistic values

    df = total_n - 2; % Degrees of freedom
    t_thresh = abs(tinv(alpha/2, df));    % Inverse Cumulative Distribution Function
   

    % 1. Calculate Observed Statistics
    obs_t_map = calculate_t_map(maps1, maps2);
    obs_stats = get_map_metrics(obs_t_map, t_thresh);

    % 2. Initialize Null Distributions
    null_dist = struct('pixel', zeros(n_perm, 1), 'size', zeros(n_perm, 1), 'mass', zeros(n_perm, 1));

    fprintf('Running %d permutations (df=14, t_thresh=%.2f)...\n', n_perm, t_thresh);
    tic;
    for i = 1:n_perm
        shuffled_idx = randperm(total_n);
        perm_t_map = calculate_t_map(all_maps(:,:,shuffled_idx(1:n1)), ...
                                     all_maps(:,:,shuffled_idx(n1+1:end)));
        m = get_map_metrics(perm_t_map, t_thresh);
        
        null_dist.pixel(i) = m.max_pixel;
        null_dist.size(i)  = m.max_size;
        null_dist.mass(i)  = m.max_mass;
        
        if mod(i, 500) == 0, fprintf('Iteration %d/%d (Elapsed: %.1fs)\n', i, n_perm, toc); end
    end

    % 3. Calculate P-values
    p_vals.pixel = sum(null_dist.pixel >= obs_stats.max_pixel) / n_perm;
    p_vals.size  = sum(null_dist.size >= obs_stats.max_size) / n_perm;
    p_vals.mass  = sum(null_dist.mass >= obs_stats.max_mass) / n_perm;

    % 4. Visualization
    plot_results(null_dist, obs_stats, p_vals, obs_t_map, t_thresh, Suc);
end

% --- Internal Helper Functions ---

function t_map = calculate_t_map(g1, g2)
    m1 = mean(g1, 3); m2 = mean(g2, 3);
    v1 = var(g1, 0, 3); v2 = var(g2, 0, 3);
    n1 = size(g1, 3); n2 = size(g2, 3);
    t_map = (m2 - m1) ./ sqrt(v1/n1 + v2/n2); % Positive t = Enriched Higher
    t_map(isnan(t_map)) = 0;
end


function stats = get_map_metrics(t_map, t_thresh)
    % Physical Scaling: 4mm / 201 pixels
    mm_per_pixel_side = 4 / 201; 
    mm2_per_pixel = mm_per_pixel_side^2; % approx 0.000396
    
    binary_map = abs(t_map) > t_thresh;
    stats.max_pixel = max(abs(t_map(:))); % T-stat is dimensionless
    
    CC = bwconncomp(binary_map);
    if CC.NumObjects > 0
        pixel_counts = cellfun(@numel, CC.PixelIdxList);
        
        % Convert Size to mm^2
        stats.max_size = max(pixel_counts) * mm2_per_pixel;
        
        % Convert Mass to mm^2 (Weighted by T-intensity)
        cluster_masses = cellfun(@(idx) sum(abs(t_map(idx))), CC.PixelIdxList);
        stats.max_mass = max(cluster_masses) * mm2_per_pixel;
    else
        [stats.max_size, stats.max_mass] = deal(0, 0);
    end
end

function plot_results(null_dist, obs, p, t_map, t_thresh, Suc)
    figure('Color', 'w', 'Position', [210 500 400 400], 'Name', 'Figure 2C');
    
    % 2D Statistical Parametric Map (SPM)

    plDIM = [4 4];
    imDIM = [201 201];
    pRES = plDIM./(imDIM-1);
    Xref = -plDIM(1):pRES(1):0;
    Yref = -plDIM(2)/2:pRES(2):plDIM(2)/2;
    
    % Cmap = fire(256);
    Cmap = parula(256);
 
    winStep = pRES(1);

     % Symmetrical limits centered at zero
    max_val = max(abs(t_map(:)));
    if max_val < t_thresh, max_val = t_thresh + 1; end % Ensure scale shows threshold

    plot3(Suc.xySensory(:,1)+winStep/2, Suc.xySensory(:,2)+winStep/2, 3*ones(1,size(Suc.xySensory,1)), 'w--')
    hold on
    plot3(Suc.xyMotor(:,1)+winStep/2, Suc.xyMotor(:,2)+winStep/2, 3*ones(1,size(Suc.xyMotor,1)), 'w--')

    surf(Xref, Yref, t_map);
    shading flat

    % Add Contours to indicate the Threshold of Significance
    % This draws a black line around any region where |t| > t_thresh
    % 1. Define your minimum cluster size in mm^2 (e.g., 0.05 mm^2)
    min_cluster_area_mm2 = 0.01; 
    
    % 2. Convert that threshold to pixels
    mm2_per_pixel = (4/201)^2;
    min_cluster_size_pixels = min_cluster_area_mm2 / mm2_per_pixel;
    
    % 3. Create a binary mask of pixels exceeding the T-threshold
    binary_map = abs(t_map) > t_thresh;
    
    % 4. Filter the binary mask: Keep only clusters larger than the pixel threshold
    % 'bwareaopen' removes all connected components (objects) that have fewer than P pixels
    clean_binary_map = bwareaopen(binary_map, round(min_cluster_size_pixels));
    
    % 5. Apply this mask to the t-map so only large clusters remain
    % We multiply the t_map by the logical mask (0s and 1s)
    filtered_t_map = abs(t_map) .* clean_binary_map;

    [~, hContour] = contour3(Xref, Yref, filtered_t_map, [t_thresh t_thresh], 'Color', 'k', 'LineWidth', 1.5);

    

    hold off
    set(gca, 'XLim', [-3.5 Xref(end)]+winStep/2, 'YLim', [Yref(1) Yref(end)]+winStep/2, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'XTick', (Xref(1):0.5:Xref(end))+winStep/2, ...
        'YTick', (Yref(1):0.5:Yref(end))+winStep/2, 'XTickLabel',  Xref(1):0.5:Xref(end),  'YTickLabel', Yref(1):0.5:Yref(end), 'XGrid', 'on', 'YGrid', 'on', 'Box', 'off', 'CLim', [-max_val, max_val])
    xlabel('Medio-lateral position (mm)')
    ylabel('Antero-posterior position (mm)')
    view(0,90)
    axis square
    colormap(gca, Cmap);

    % Formatting the Colorbar
    hCB = colorbar;
    hCB.Label.String = 'T-Value (Signal-to-Noise Ratio)';

    title(['T-Map: EE vs CH', newline, 'Black Contour: |t| > ', num2str(t_thresh)]);

%     How to interpret the scale:
% 
%     A value of 0: There is no difference between the Control and Enriched means.
% 
%     A value of +2.14: The Enriched group is 2.14 standard errors higher than the Control group. At sample size N=16, this is the "threshold of significance" (p=0.05).
% 
%     A value of -2.14: The Control group is 2.14 standard errors higher than the Enriched group.


     % Plot histograms
    figure('Color', 'w', 'Position', [10 300 800 400], 'Name', 'Figure 2D');
    fields = {'size', 'mass'};
    labels = {'Max Cluster Size', 'Max Cluster Mass'};

   
    subplot(121)
    h = histogram(null_dist.size, 0:0.05:1, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'w');
    h.Normalization = 'probability';
    hold on;
    line([obs.('max_size') obs.('max_size')], ylim, ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    xlabel('Maximum cluster size (mm^2)');
    text(0.5, 0.2, sprintf('p = %.3f', p.('size')))
    set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'off')

    subplot(122)
    h = histogram(null_dist.mass, 0:0.2:4, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'w');
    h.Normalization = 'probability';
    hold on;
    line([obs.('max_mass') obs.('max_mass')], ylim, ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    xlabel('Maximum cluster mass (mm^{2}.t)');
    text(1.5, 0.2, sprintf('p = %.3f', p.('mass')))
    set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'off')
    


end