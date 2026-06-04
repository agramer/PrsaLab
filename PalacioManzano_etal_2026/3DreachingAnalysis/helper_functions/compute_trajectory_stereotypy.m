function [pca_stereotypy, xyzAll_common] = compute_trajectory_stereotypy(XYZ, Ncomp)
    % Computes PCA-based stereotypy
    %
    % Inputs arguments:
    %   XYZ - 1x3 cell array containing X, Y, Z coordinate matrices (Ntrials x Nsamples each)
    %         
    %
    % Outputs arguments:
    %   pca_stereotypy - PCA-based stereotypy score (% variance in first Ncomp components)
        
     
        
    % Align trajectories for PCA
    aligned_trajs = align_trajectories(XYZ);
    
    % Concatenate all dimensions
    data_matrix = [aligned_trajs{1}, aligned_trajs{2}, aligned_trajs{3}];
    
    % Remove any remaining NaN rows
    valid_rows = ~any(isnan(data_matrix), 2);
    data_matrix = data_matrix(valid_rows, :);
    
    % Perform PCA
    [coeff, score, latent, ~, explained] = pca(data_matrix);
    
    % PCA stereotypy as cumulative variance in first few components
    n_components = min(Ncomp, length(explained));
    pca_stereotypy = sum(explained(1:n_components));
    
    xyzAll_common = aligned_trajs';
   
end

function aligned_trajs = align_trajectories(XYZ)
    % Align trajectories to common length using interpolation
 
    n_trials = size(XYZ{1}, 1);
    common_length = 50;
    aligned_trajs = cell(3, 1);
    for dim = 1:3
        aligned_trajs{dim} = nan(n_trials, common_length);
    end
    
    % Interpolate each trial to common length
    for trial = 1:n_trials
        for dim = 1:3
            valid_data = XYZ{dim}(trial, ~isnan(XYZ{dim}(trial, :)));
            if length(valid_data) > 1  
                aligned_trajs{dim}(trial, :) = interp1(1:length(valid_data), ...
                    valid_data, linspace(1, length(valid_data), common_length));
            end
        end
    end
end


