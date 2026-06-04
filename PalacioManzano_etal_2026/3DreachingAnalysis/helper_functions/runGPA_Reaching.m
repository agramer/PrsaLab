function xyzAll_aligned = runGPA_Reaching(xyzAll, Nrep)
    % GPA_TRAJECTORIES  Generalized Procrustes Analysis for 3D reach trajectories.
    %
    % Inputs:
    %   xyzAll      : 1x3 cell array containing X, Y, Z coordinate matrices (Ntrials x Nsamples each)
    %   Nrep        : number of iterations to get the mean template
    %
    % Output:
    %   xyzAll_aligned : 1x3 cell array in same format as xyzAll, GPA-aligned
    
    Ntrials  = size(xyzAll{1}, 1);
    Nsamples = size(xyzAll{1}, 2);
    
    trajs = zeros(Ntrials, Nsamples, 3);
    for d = 1:3
        trajs(:, :, d) = xyzAll{d};
    end
    
    % Remove translation
    for i = 1:Ntrials
        for d = 1:3
            trajs(i, :, d) = trajs(i, :, d) - mean(trajs(i, :, d));
        end
    end
    
    % Find most typical trial (medoid) first
    [medoid_trial, ~] = find_medoid_trial(trajs);
    
    % Initialize template as the medoid trajectory 
    template = squeeze(trajs(medoid_trial, :, :));   
    
    % Run GPA
    for rep = 1:Nrep
    
        % Align each trajectory to the current template
        for i = 1:Ntrials
            T = squeeze(trajs(i, :, :));   
            [~, T_aligned] = procrustes(template, T, 'Scaling', 0, 'Reflection', 0);
            trajs(i, :, :) = T_aligned;
        end
    
        % Update template as the mean of all aligned trajectories
        template = squeeze(mean(trajs, 1));   % (Nsamples x 3)
    
        % Re-center the template
        for d = 1:3
            template(:, d) = template(:, d) - mean(template(:, d));
        end
    
    end
    
    % Normalize 
    trajs = trajs - mean(trajs(:,1,:), 1);
    
    xyzAll_aligned = cell(1, 3);
    for d = 1:3
        xyzAll_aligned{d} = trajs(:, :, d);
    end

end


function [medoid_trial, D] = find_medoid_trial(traj)

    Ntrials = size(traj, 1);
    
    % Compute pairwise Procrustes distances
    D = zeros(Ntrials, Ntrials);
    for i = 1:Ntrials
        for j = i+1:Ntrials
            T_i = squeeze(traj(i, :, :));   
            T_j = squeeze(traj(j, :, :));   
            d = procrustes(T_i, T_j, 'Scaling', 0, 'Reflection', 0);
            D(i, j) = d;
            D(j, i) = d;                   
        end
    end
    
    % Medoid is at minimal distance to all other trials
    total_distances  = sum(D, 2);          
    [~, medoid_trial] = min(total_distances);

end