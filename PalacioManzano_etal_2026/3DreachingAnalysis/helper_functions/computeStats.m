function [ReachStats, GraspStats, ReturnStats, FullTrajStats, xyzReach, xyzGrasp, xyzReturn, mvtVars ] = computeStats(xyzAll, velAll, mvtVars)

Ntrials = size(xyzAll{1,1},1);
ReachStats = zeros(Ntrials, 4);     % distance, peak velocity, duration (in samples), path length
GraspStats = zeros(Ntrials, 4);     % distance, peak velocity, duration (in samples), path length
ReturnStats = zeros(Ntrials, 4);    % distance, peak velocity, duration (in samples), path length
FullTrajStats = zeros(Ntrials, 2);  % path length, peak velocity

xyzReach = cell(size(xyzAll));
xyzGrasp = cell(size(xyzAll));
xyzReturn = cell(size(xyzAll));
xyzFullTraj = cell(size(xyzAll));

for dim=1:3
    maxReachSamps = max(mvtVars(:,2)-mvtVars(:,1))+1;
    xyzReach{1,dim} = nan(Ntrials,maxReachSamps);

    maxGraspSamps = max(mvtVars(:,4)-mvtVars(:,3))+1;
    xyzGrasp{1,dim} = nan(Ntrials,maxGraspSamps);

    maxReturnSamps = max(mvtVars(:,6)-mvtVars(:,5))+1;
    xyzReturn{1,dim} = nan(Ntrials,maxReturnSamps);   
end

for tr=1:Ntrials

    NsampsReach = mvtVars(tr,2)-mvtVars(tr,1)+1;
    NsampsGrasp = mvtVars(tr,4)-mvtVars(tr,3)+1;
    NsampsReturn = mvtVars(tr,6)-mvtVars(tr,5)+1;
    NsampsFullTraj = mvtVars(tr,6)-mvtVars(tr,1)+1;
    velReach = zeros(3,NsampsReach);
    velGrasp = zeros(3,NsampsGrasp);
    velReturn = zeros(3,NsampsReturn);
    velFullTraj = zeros(3,NsampsFullTraj);
    
    
    for dim=1:3
        xyzReach{1,dim}(tr,1:NsampsReach) = xyzAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,2));
        xyzGrasp{1,dim}(tr,1:NsampsGrasp) = xyzAll{1,dim}(tr,mvtVars(tr,3):mvtVars(tr,4));
        xyzReturn{1,dim}(tr,1:NsampsReturn) = xyzAll{1,dim}(tr,mvtVars(tr,5):mvtVars(tr,6));
        xyzFullTraj{1,dim}(tr,1:NsampsFullTraj) = xyzAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,6));

        velReach(dim, :) =  velAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,2));
        velGrasp(dim, :) =  velAll{1,dim}(tr,mvtVars(tr,3):mvtVars(tr,4));
        velReturn(dim, :) =  velAll{1,dim}(tr,mvtVars(tr,5):mvtVars(tr,6));
        velFullTraj(dim, :) =  velAll{1,dim}(tr,mvtVars(tr,1):mvtVars(tr,6));
    end
    ReachStats(tr, 3) = NsampsReach;        % duration
    GraspStats(tr, 3) = NsampsGrasp;        % duration
    ReturnStats(tr, 3) = NsampsReturn;      % duration
    FullTrajStats(tr, 3) = NsampsFullTraj;  % duration
    
    xyzFullStart = [xyzFullTraj{1,1}(tr,1);xyzFullTraj{1,2}(tr,1);xyzFullTraj{1,3}(tr,1)];
    xyzFullEnd = [xyzFullTraj{1,1}(tr,NsampsFullTraj);xyzFullTraj{1,2}(tr,NsampsFullTraj);xyzFullTraj{1,3}(tr,NsampsFullTraj)];
    FullTrajStats(tr, 1) = vecnorm(xyzFullEnd-xyzFullStart);     % distance
    [FullTrajStats(tr, 4), ~] = getPathLength(xyzFullTraj,tr);    % path length
    [FullTrajStats(tr, 2), ~] = max(sqrt(velFullTraj(1,:).^2+velFullTraj(2,:).^2+velFullTraj(3,:).^2));  % peak velocity

    
    xyzReachStart = [xyzReach{1,1}(tr,1);xyzReach{1,2}(tr,1);xyzReach{1,3}(tr,1)];
    xyzReachEnd = [xyzReach{1,1}(tr,NsampsReach);xyzReach{1,2}(tr,NsampsReach);xyzReach{1,3}(tr,NsampsReach)];
    ReachStats(tr, 1) = vecnorm(xyzReachEnd-xyzReachStart);     % distance
    [ReachStats(tr, 2), mvtVars(tr,7)] = max(sqrt(velReach(1,:).^2+velReach(2,:).^2+velReach(3,:).^2));  % peak velocity
    mvtVars(tr,7) = mvtVars(tr,7) + mvtVars(tr,1) - 1;      % peak velocity time
    [ReachStats(tr, 4), midInd] = getPathLength(xyzReach,tr);    % path length

    
    xyzGraspStart = [xyzGrasp{1,1}(tr,1);xyzGrasp{1,2}(tr,1);xyzGrasp{1,3}(tr,1)];
    xyzGraspEnd = [xyzGrasp{1,1}(tr,NsampsGrasp);xyzGrasp{1,2}(tr,NsampsGrasp);xyzGrasp{1,3}(tr,NsampsGrasp)];
    GraspStats(tr, 1) = vecnorm(xyzGraspEnd-xyzGraspStart);     % distance
    [GraspStats(tr, 2), mvtVars(tr,8)] = max(sqrt(velGrasp(1,:).^2+velGrasp(2,:).^2+velGrasp(3,:).^2));  % peak velocity
    mvtVars(tr,8) = mvtVars(tr,8) + mvtVars(tr,3) - 1;      % peak velocity time
    [GraspStats(tr, 4), midInd] = getPathLength(xyzGrasp,tr);    % path length
    
    xyzReturnStart = [xyzReturn{1,1}(tr,1);xyzReturn{1,2}(tr,1);xyzReturn{1,3}(tr,1)];
    xyzReturnEnd = [xyzReturn{1,1}(tr,NsampsReturn);xyzReturn{1,2}(tr,NsampsReturn);xyzReturn{1,3}(tr,NsampsReturn)];
    ReturnStats(tr, 1) = vecnorm(xyzReturnEnd-xyzReturnStart);      % distance
    [ReturnStats(tr, 2), mvtVars(tr,9)] = max(sqrt(velReturn(1,:).^2+velReturn(2,:).^2+velReturn(3,:).^2));  % peak velocity
    mvtVars(tr,9) = mvtVars(tr,9) + mvtVars(tr,5) - 1;      % peak velocity time
    [ReturnStats(tr, 4), midInd] = getPathLength(xyzReturn, tr);    % path length

end