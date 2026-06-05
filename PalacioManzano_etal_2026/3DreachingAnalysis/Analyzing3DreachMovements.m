


fpath = 'm3061_Session1_PRE.mat';   % insert name of the file from raw_data folder to process
S= load(fpath);





addpath(genpath('helper_functions'))
addpath('raw_data')

velThrVxVy = [20 20];
xStartStop = -14;
normFlag = 1;

Coords3D = S.Coords3D;
JointCoords = S.JointCoords;
datTimings = S.datTimings;
clear S
percCorr = 100*sum(datTimings(:,10), 'omitnan')/nnz(~isnan(datTimings(:,10)));



%%%%%%%%%%%%%%%%%%%%%%%
% Pre-process data
%%%%%%%%%%%%%%%%%%%%%%%

Ntrials = size(JointCoords,1);
% Exclude invalid trials
for tr=1:Ntrials
    if ~isempty(JointCoords{tr,1}) && JointCoords{tr,1}(2,1,1)==20
        JointCoords{tr,1}=[];
        Coords3D{tr,1}=[];
    end
end
indTr = find(~cellfun(@isempty, Coords3D));
indTr = indTr';
Ntrials = size(indTr,2);

% load stereo camera calibration transform
S=load('250530_CalibCoords.mat');
tform = S.tform;
clear S;

maxSamps = max(cellfun(@(x) size(x,3), Coords3D));
Fs = 100;   % Sampling rate (Hz, fps)
periTime = 0:1/Fs:(maxSamps-1)/Fs;

xyzAll = cell(1,3);
for dim=1:3
    xyzAll{1,dim} = zeros(Ntrials,numel(periTime));
end


for cnt = 1:Ntrials 
    trInd = indTr(cnt);

    xyz = nan(3, maxSamps);
    Nsamps = size(Coords3D{trInd,1},3);
    xyz(:,1:Nsamps) = reshape(Coords3D{trInd,1}, [3 Nsamps]);
    xyz=xyz';
    
    % Low pass filter the data to correct for tracking noise
    ind = find(isnan(xyz(:,1)),1);
    if isempty(ind)
        xyz = sgolayfilt(xyz, 2, 9);
    else
        xyz(1:ind-1,:) = sgolayfilt(xyz(1:ind-1,:), 2, 9);
    end

    xyz = transformPointsForward(tform,xyz);
    xyz = xyz';

    for dim=1:3
        xyzAll{1,dim}(cnt,:) = xyz(dim,:);
    end
end

velAll = cell(size(xyzAll));
for dim=1:3
    velAll{1,dim} = gradient(xyzAll{1,dim}, 1/Fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Reach, Graps, Retract components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntrials = size(xyzAll{1,1},1);
Nsamps = size(xyzAll{1,1},2);
midSampRatio = 3/4;         % look for maximal extension within the first 3/4 of samples (to avoid detecting second reaches late in the trial)
midSampLim = round(Nsamps*midSampRatio);

mvtVars = nan(Ntrials,10);       % reach onset, reach end, grasp onset, grasp end, return onset, return end, reach midpoint time, grasp midpoint time, return midpoint time, trial index  
mvtVars(:,10) = 1:Ntrials;       % to keep track of outliers


                                  
maxSamps = 0;
for tr=1:Ntrials
    [~,indMaxX] = max(xyzAll{1,1}(tr,1:midSampLim));
    [pks, locs] = findpeaks(velAll{1,1}(tr, indMaxX:-1:1));

    if ~isempty(pks)
       [~,ind] = max(pks);
       indPeakVx = indMaxX - locs(ind) + 1;
    else
       indPeakVx = indMaxX;
    end
    indStart = find(velAll{1,1}(tr, indPeakVx:-1:1)<=0, 1);
    if isempty(indStart)
        [pks, locs] = findpeaks(-velAll{1,1}(tr, indPeakVx:-1:1));  % find throughs
        if ~isempty(pks)
            [~,ind] = max(pks);
            indStart = indPeakVx - locs(ind) + 1;
        else
            indStart = 1;
        end
    else
        indStart = indPeakVx - indStart + 1;
    end
    indStop = find(xyzAll{1,1}(tr,indMaxX:end)<xStartStop, 1);
    indStop = indMaxX + indStop - 1;
    if isempty(indStop)
        indStop = find(velAll{1,1}(tr,Nsamps:-1:indMaxX)<0, 1);
        indStop = Nsamps - indStop + 1;
    end
    if isempty(indStop) || indStop>Nsamps
        indStop = Nsamps;
    end
    

    indReachStart = find(velAll{1,1}(tr, indStart:indStop)>velThrVxVy(1), 1);
    if isempty(indReachStart)
        indReachStart = indStart;
    else
        indReachStart = indStart + indReachStart - 1;
    end
    mvtVars(tr,1) =  indReachStart; % start of reach phase
 

    [~,indTroughVx] = min(velAll{1,1}(tr, indMaxX:indStop));
    indTroughVx = indMaxX + indTroughVx - 1;

    [~, locs] = findpeaks(-velAll{1,2}(tr, indPeakVx:indTroughVx));  % find throughs
    locsT = locs + indPeakVx - 1;
    [~, ind] = max(xyzAll{1,1}(tr,locsT)); % find peak Vy with most extreme X position
    indPeakVy = locs(ind);
    
    indPeakVy = indPeakVx + indPeakVy - 1;
    indGraspStart = find(velAll{1,2}(tr, indPeakVy:-1:indReachStart)>-velThrVxVy(2), 1);
    if isempty(indGraspStart)
        indGraspStart = indPeakVx;
    else
        indGraspStart = indPeakVy - indGraspStart + 1;
    end
    mvtVars(tr,3) =  indGraspStart;  % start of grasp phase 


    indGraspEnd = find(velAll{1,2}(tr, indPeakVy:indStop)>-velThrVxVy(2), 1);
    if isempty(indGraspEnd)
        indGraspEnd = indTroughVx;
    else
        indGraspEnd = indPeakVy + indGraspEnd - 1;
    end
    mvtVars(tr,4) =  indGraspEnd;  % end of grasp phase


    indStop = find(velAll{1,1}(tr, indTroughVx:end)>=0, 1);
    if isempty(indStop)
        indStop = Nsamps;
    else
        indStop = indTroughVx + indStop - 1;
    end

    indReturnEnd = find(velAll{1,1}(tr, indStop:-1:indStart)<-velThrVxVy(1), 1);
    if isempty(indReturnEnd)
        indReturnEnd = indStop;
    else
        indReturnEnd = indStop - indReturnEnd + 1;
    end
    mvtVars(tr,6) =  indReturnEnd; % end of return phase


    for dim=1:3
        dat = xyzAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,6));
        xyzAll{1,dim}(tr,:) = nan;
        xyzAll{1,dim}(tr,1:mvtVars(tr,6)-mvtVars(tr,1)+1) = dat;

        dat = velAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,6));
        velAll{1,dim}(tr,:) = nan;
        velAll{1,dim}(tr,1:mvtVars(tr,6)-mvtVars(tr,1)+1) = dat;
    end

    maxSamps = max(maxSamps, mvtVars(tr,6)-mvtVars(tr,1)+1);
    mvtVars(tr,1:end-1) = mvtVars(tr,1:end-1) - mvtVars(tr,1) + 1;

end

mvtVars(:,2) = mvtVars(:,3);
mvtVars(:,5) = mvtVars(:,4);
mvtVars = round(mvtVars);

for dim=1:3
    xyzAll{1,dim}(:,maxSamps+1:end)=[];
    velAll{1,dim}(:,maxSamps+1:end)=[];
end  

% Standardize starting positions
for dim=1:3
    MvtStats{1,2}(:,dim) = xyzAll{1,dim}(:,1);  % starting xyz coordinates of reaching
    xyzAll{1,dim} = xyzAll{1,dim} - mean(xyzAll{1,dim}(:,1), 'omitnan');
end


% Remove outliers based on path lenghts
pathLengths = zeros(Ntrials,3);
for tr=1:Ntrials
    xyzReach = cellfun(@(x) x(:,mvtVars(tr,1):mvtVars(tr,2)), xyzAll, 'UniformOutput', false);
    pathLengths(tr,1) = getPathLength(xyzReach, tr);
    xyzGrasp = cellfun(@(x) x(:,mvtVars(tr,3):mvtVars(tr,4)), xyzAll, 'UniformOutput', false);
    pathLengths(tr,2) = getPathLength(xyzGrasp, tr);
    xyzReturn = cellfun(@(x) x(:,mvtVars(tr,5):mvtVars(tr,6)), xyzAll, 'UniformOutput', false);
    pathLengths(tr,3) = getPathLength(xyzReturn, tr);
end
indOut =  isoutlier(pathLengths(:,1), 'gesd') | isoutlier(pathLengths(:,3), 'gesd');
disp(['Number of outliers: ' num2str(sum(indOut))]);
xyzAll = cellfun(@(x) x(~indOut, :), xyzAll, 'UniformOutput', false);
velAll = cellfun(@(x) x(~indOut, :), velAll, 'UniformOutput', false);
mvtVars(indOut,:) = [];

% Re-standardize 
for dim=1:3
    xyzAll{1,dim} = xyzAll{1,dim} - mean(xyzAll{1,dim}(:,1), 'omitnan');
end

[ReachStats,GraspStats, ReturnStats, ~, ~, ~, ~, mvtVars] = computeStats(xyzAll, velAll, mvtVars);

ampThr = 2; 
indOut = isoutlier(ReachStats(:,1), 'gesd') | ReachStats(:,1)<ampThr | isoutlier(GraspStats(:,1), 'gesd') | GraspStats(:,1)<ampThr | isoutlier(ReturnStats(:,1), 'gesd') | ReturnStats(:,1)<ampThr;
xyzAll = cellfun(@(x) x(~indOut, :), xyzAll, 'UniformOutput', false);
velAll = cellfun(@(x) x(~indOut, :), velAll, 'UniformOutput', false);
mvtVars(indOut,:) = [];
MvtStats{1,2}(indOut,:) = [];
for dim=1:3
    xyzAll{1,dim} = xyzAll{1,dim} - mean(xyzAll{1,dim}(:,1), 'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute kinematics of each component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ReachStats, GraspStats, ReturnStats, FullTrajStats, xyzReach, xyzGrasp, xyzReturn, mvtVars] = computeStats(xyzAll, velAll, mvtVars);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PCA-based stereotypy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ncomp = 3;              % number of PCA components for stereotypy calculation
[pca_stereotypy, xyzAll_common] = compute_trajectory_stereotypy(xyzAll, Ncomp);

xyzAll_resampled = arc_length_reparam(xyzAll_common, 50);       % arc-length re-parameterization
xyzAll_aligned = runGPA_Reaching(xyzAll_resampled, 100);     % generalized procrustes transformation (translation + rotation)

[pca_stereotypy_aligned, ~] = compute_trajectory_stereotypy(xyzAll_aligned, Ncomp);

disp(['Stereotypy index = ' num2str(pca_stereotypy)])
disp(['Stereotypy index (post-correction) = ' num2str(pca_stereotypy_aligned)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot kinematics and get main-sequence slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varNames = {'Movement amplitude (mm)'; 'Peak velocity (mm/s)'; 'Movement duration (samples)'; 'Path length (mm)'};

cols = {[[0.8 0.8 0.8]; [0.45 0.7 0.9]; [0 0.32 0.8]];[[0.8 0.8 0.8]; [0.9 0.5 0.5]; [0.9 0.2 0.2]];[[0.8 0.8 0.8]; [0.7 0.85 0.55]; [0.4 0.6 0.2]]};

% col = [1 0.7 0.5];  % orange
col = [0.9 0.9 0.9]; % grey

regStats = zeros(4, 6);


trialInds = [];


varIDa = 1;     % 1: distance, 2: peak velocity, 3: duration (in samples), 4: path length
varIDb = 2;

figure('Color', 'w')
hAx=gca;
plotTrajectories(xyzAll, hAx, [1 0.7 0.5], normFlag, trialInds)  

figure('Color', 'w')
subplot(231)
hAx=gca;
% plot the two separately to avoid artefacts when exporting figure to ps
plotTrajectories(xyzAll, hAx, col, normFlag, trialInds)  
plotTrajectories(xyzReach, hAx, cols{1,1}(3,:), normFlag, trialInds)
subplot(234)
hAx=gca;
regStats(1,:) = plotTrialKinematics(hAx, ReachStats,  cols{1,1}(3,:), varNames, varIDa, varIDb);
subplot(232)
hAx=gca;
plotTrajectories(xyzAll, hAx, col, normFlag, trialInds)
plotTrajectories(xyzGrasp, hAx,  cols{2,1}(3,:), normFlag, trialInds)
subplot(235)
hAx=gca;
regStats(2,:) = plotTrialKinematics(hAx, GraspStats,  cols{2,1}(3,:),  varNames, varIDa, varIDb);
subplot(233)
hAx=gca;
plotTrajectories(xyzAll, hAx, col, normFlag, trialInds)
plotTrajectories(xyzReturn, hAx,  cols{3,1}(3,:), normFlag, trialInds)
subplot(236)
hAx=gca;
regStats(3,:) = plotTrialKinematics(hAx, ReturnStats,  cols{3,1}(3,:), varNames, varIDa, varIDb);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot elevation and azimuth angles and compute means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzimuthMeanAngles = zeros(1,3);
ElevationhMeanAngles = zeros(1,3);

[AzimuthMeanAngles(1,1), ElevationhMeanAngles(1,1)] = computeAngles(xyzReach, 1);
[AzimuthMeanAngles(1,2), ElevationhMeanAngles(1,2)] = computeAngles(xyzGrasp, 2);
[AzimuthMeanAngles(1,3), ElevationhMeanAngles(1,3)] = computeAngles(xyzReturn, 3);

f1=figure('Color','w', 'Name', 'Elevation angle');
f2=figure('Color','w', 'Name', 'Azimuth angle');
for k=1:3
    figure(f1);
    subplot(1,3,k)
    hold on
    plotArc(-90, 90, 1, gca)
    axis equal
    set(gca, 'XLim', [-0.1 1.1], 'YLim', [-1.1 1.1], 'Visible', 'off')

    figure(f2);
    subplot(1,3,k)
    hold on
    plotArc(-180, 180, 1, gca)
    axis equal
    set(gca, 'XLim', [-1.1 1.1], 'YLim', [-1.1 1.1],  'Visible', 'off')
end

for comp=1:3
    figure(f1);
    subplot(1,3,comp)
    theta = ElevationhMeanAngles(1,comp);
    x = [0, cosd(theta)];    % cosd, sind use degrees
    y = [0, sind(theta)];
    plot(x, y, 'Color', cols{comp,1}(2,:), 'LineWidth', 2);
    
    figure(f2);
    subplot(1,3,comp)
    theta = -AzimuthMeanAngles(1,comp)+90;
    x = [0, cosd(theta)];    % cosd, sind use degrees
    y = [0, sind(theta)];
    plot(x, y, 'Color', cols{comp,1}(2,:), 'LineWidth', 2);
end

disp(['Reach elevation=' num2str(ElevationhMeanAngles(1,1)), '  Reach azimuth=' num2str(AzimuthMeanAngles(1,1))])
disp(['Grasp elevation=' num2str(ElevationhMeanAngles(1,2)), '  Grasp azimuth=' num2str(AzimuthMeanAngles(1,2))])
disp(['Return elevation=' num2str(ElevationhMeanAngles(1,3)), '  Return azimuth=' num2str(AzimuthMeanAngles(1,3))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample-wise dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = xyzAll_common{1}';
Y = xyzAll_common{2}';
Z = xyzAll_common{3}';

mu_X = mean(X, 2);
mu_Y = mean(Y, 2);
mu_Z = mean(Z, 2);
dist  = sqrt((X - mu_X).^2 + (Y - mu_Y).^2 + (Z - mu_Z).^2);
mean_dist = mean(dist, 2);


Nbins = 5;
Nsamps = size(mean_dist,1);

figure('Color','w')
plot(mean_dist, 'Color', 'k')

plot(mean_dist, 'Color', 'k', 'LineWidth',2)
set(gca, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'off', 'XTick', 0.5:Nsamps/Nbins:Nsamps+0.5, 'XTickLabel', 0:1/Nbins:1, 'XGrid', 'on', 'XLim', [0 Nsamps+1])
xlabel('Normalized path')
ylabel('Sample-wise dispersion')

