
addpath('helper_functions')
addpath('raw_data')

velThrVxVy = [20 20];
xStartStop = -14;
normFlag = 1;

procrustFlag = 0;

fpath = [];
if isempty(fpath)
    S=load('m3061_Session1_PRE.mat');
else
    S= load(fpath);
end

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


%%%%%%%%
% Azimuth and Elevation
- plot and output mean values

%%%%
Path efficiency

