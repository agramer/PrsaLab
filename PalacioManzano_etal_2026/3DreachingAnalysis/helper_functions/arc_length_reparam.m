function xyzAll_resampled = arc_length_reparam(xyzAll, N_out)
% Resample 3D trajectories at uniform arc-length intervals.
%
% Input arguments:
%   xyzAll  : 1x3 cell array containing X, Y, Z coordinate matrices (Ntrials x Nsamples each)
%   N_out   : number of output samples
%
% Output arguments:
%   xyzAll_resampled : 1x3 cell array in same format as xyzAll

Ntrials  = size(xyzAll{1}, 1);

xyzAll_resampled = cell(1, 3);
for d = 1:3
    xyzAll_resampled{d} = zeros(Ntrials, N_out);
end

for i = 1:Ntrials

    % Extract (Nsamples x 3) trajectory for this trial
    traj = zeros(size(xyzAll{1}, 2), 3);
    for d = 1:3
        traj(:, d) = xyzAll{d}(i, :)';
    end

    % Compute cumulative arc length
    diffs       = diff(traj, 1, 1);                    
    seg_lengths = sqrt(sum(diffs.^2, 2));              
    arc         = [0; cumsum(seg_lengths)];            
    arc         = arc / arc(end);                      % normalize to [0, 1]

    % Remove duplicate arc values to ensure valid interpolation
    [arc_unique, idx] = unique(arc);
    t_uniform = linspace(0, 1, N_out)';

    % Interpolate 
    for d = 1:3
        xyzAll_resampled{d}(i, :) = interp1(arc_unique, traj(idx, d), t_uniform, 'pchip');
    end

end

end

