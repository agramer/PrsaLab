function [PDs, tuneAct] = getMesoMapPD(ANAframes_DIR, msIDs)


Ndirs = 8;
angles = 0:pi/4:2*pi-pi/4;
angles_deg = 0:45:315;
Nms = numel(msIDs); % number of mice or sessions
maxNormFlag = 1;
tuneAct = zeros(Nms,Ndirs);

PDs = zeros(Nms,1); % preferred directions of the mesoscale activation maps

for k=1:Nms
    for dr=1:Ndirs
        actFR = ANAframes_DIR{dr,1}(:,:,msIDs(k));
        tuneAct(k,dr) = max(actFR(:));
    end
end


for k=1:Nms
    if maxNormFlag
        tuneAct(k,:) = (tuneAct(k,:) - min(tuneAct(k,:)));
        tuneAct(k,:) = tuneAct(k,:) / max(tuneAct(k,:));
    end
    [X,Y] = pol2cart(angles,tuneAct(k,:));
    X = [X X(1)];
    Y = [Y Y(1)];
    pIN = polyshape(X,Y);
    [xC,yC] = centroid(pIN);
    [th,~] = cart2pol(xC,yC);
    PDs(k,1) = mod(rad2deg(th), 360);
end


