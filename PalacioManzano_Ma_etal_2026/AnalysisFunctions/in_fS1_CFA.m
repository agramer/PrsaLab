function R = in_fS1_CFA(CenterXY, S)
% IN_FS1_CFA  Which cells lie within the fS1 and CFA boundaries defined by the
% Ueno et al. (2018) contours.


valid = all(~isnan(CenterXY), 2);
x = CenterXY(:,1);  
y = CenterXY(:,2);
[in, on] = inpolygon(x, y, S.xySensory(:,1), S.xySensory(:,2));
R.inS1 = (in | on) & valid;
[in, on] = inpolygon(x, y, S.xyMotor(:,1), S.xyMotor(:,2));
R.inCFA = (in | on) & valid;
R.inEither = R.inS1 | R.inCFA;

R.nTotal = sum(valid);
R.nS1 = sum(R.inS1);  
R.nCFA = sum(R.inCFA);  
R.nEither = sum(R.inEither);


