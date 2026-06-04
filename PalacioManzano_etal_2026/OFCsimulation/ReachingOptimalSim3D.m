% xyzPos = ReachingOptimalSim3D(lesionPerc, lesionParamIdx)
% 
% 
% Simulate translation of a single point mass in 3D based with the OFC model
% Based on: XXX
% 
% Parameter optimization is based on a generalized LQG solution.
% Based on: Todorov, E. (2005) (see kalman_lqg_xstar.m for details)
% 
% Model paramters are described and set in load_parameters.m 
% 
% Input arguments:
%     - lesionPerc: (0 to 100) simulates lesion by downscaling OFC model
%     parameters (Default: 0)
%     - lesionPramIdx: (1 to 15) selects which OFC parameter to downscale
%     according to below definitions (Default: 8)
% 
%      1:'A', 2:'A_hat', 3:'B', 4:'B_hat', 5:'H', 6:'H_hat', 7:'L', 8:'K'
% 
% Output argument:
%     - xyzPos: 1 by 3 cell of x, y, z simulated positions (each has Ntrials x Nsamples data points) 



function xyzPos = ReachingOptimalSim3D(lesionPerc, lesionParamIdx)

addpath('helper_functions')

if nargin<1 
    lesionPerc =  0;
    lesionParamIdx = 8;
elseif nargin<2
    lesionParamIdx = 8;
end

lesionX = 1-lesionPerc/100;  
paramMultiplier = ones(1,8);
paramMultiplier(1,lesionParamIdx) = lesionX;


figure('Color', 'w', 'Position', [100 100 560 420])

load_parameters;

N = round(T_end/dt);
col = [0.7 0.7 0.7];
szX = 6*dim;
Wt = repelem(Wt,dim);
szY = size(Wt,1);


N_at_A = round(T_A/dt);         % sample number at target A
N_at_B = round((T_A+T_B)/dt);   % sample number at target B
N_at_C = N;                     % sample number at target C

A = zeros(szX,szX); 

R = eye(dim,dim)*r/(N-1);       % effort penalty matrix

D(:,:,1) = [sigma_s*diag(Wt) zeros(szY,szX-szY)];   % state dependent feedback noise

szD = size(D,3);
H = [eye(szY,szY) zeros(szY,szX-szY)];

S1  = ones(1,szX)*S1scale;                                          % additive starting point noise
E0  = ones(1,szX)*E0scale;                                          % additive internal noise
D0  = ones(1,szY)*D0scale;                                          % additive sensory noise
C0 = [zeros(1,2*dim) ones(1,dim) zeros(1,3*dim)]*C0scale;           % additive system noise 


S1 = diag(S1);
C0 = diag(C0);
E0 = diag(E0);
D0 = diag(D0);


% Three dimensional implementation
% System state: Xt = [posX; posY; posZ; velX; velY; velZ; Fx; Fy; Fz; trgA_X; trgA_Y; trgA_Z; trgB_X; trgB_Y; trgB_Z; trgC_X; trgC_Y; trgC_Z];
X1 = [p_zero(1); p_zero(2); p_zero(3); 0; 0; 0; 0; 0; 0; p_trgA(1); p_trgA(2); p_trgA(3); p_trgB(1); p_trgB(2); p_trgB(3); p_trgC(1); p_trgC(2); p_trgC(3)];     % initial state

% Vectors for computing cost matrices

% Positions at target A
px_A = [1;0;0;0;0;0;0;0;0;-1;0;0;0;0;0;0;0;0];        % Xt'*(p*p')*Xt = (position(Ttrgt)-p_start)^2;
py_A = [0;1;0;0;0;0;0;0;0;0;-1;0;0;0;0;0;0;0];
pz_A = [0;0;1;0;0;0;0;0;0;0;0;-1;0;0;0;0;0;0];

% Positions at target B
px_B = [1;0;0;0;0;0;0;0;0;0;0;0;-1;0;0;0;0;0];        % Xt'*(p*p')*Xt = (position(Ttrgt)-p_start)^2;
py_B = [0;1;0;0;0;0;0;0;0;0;0;0;0;-1;0;0;0;0];
pz_B = [0;0;1;0;0;0;0;0;0;0;0;0;0;0;-1;0;0;0];

% Positions at target C
px_C = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-1;0;0];        % Xt'*(p*p')*Xt = (position(Ttrgt)-p_start)^2;
py_C = [0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-1;0];
pz_C = [0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;-1];


vx = [0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0];                % velocity at 0
vy = [0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0];                % velocity at 0
vz = [0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0];                % velocity at 0
fx = [0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0];                % force at 0
fy = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0];                % force at 0
fz = [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0];                % force at 0


% U is 3x1 [ux;uy;uz]

szU = 3;

A(1,1) = 1;     % posX vs. posX
A(1,4) = dt;    % posX vs. velX
A(2,2) = 1;     % posY vs. posY
A(2,5) = dt;    % posY vs. velY
A(3,3) = 1;     % posZ vs. posZ
A(3,6) = dt;    % posZ vs. velZ

A(4,4) = 1-dt*b/m;  % velX vs. velX
A(4,7) = dt/m;      % velX vs. Fx
A(5,5) = 1-dt*b/m;  % velY vs. velY
A(5,8) = dt/m;      % velY vs. Fy
A(6,6) = 1-dt*b/m;  % velZ vs. velZ
A(6,9) = dt/m;      % velZ vs. Fz

A(7,7) = 1-dt/theta;    % Fx vs. Fx
A(8,8) = 1-dt/theta;    % Fy vs. Fy
A(9,9) = 1-dt/theta;    % Fz vs. Fz

for k=10:18
    A(k,k) = 1;
end


B = zeros(szX,szU);
B(:,1) = [0; 0; 0; 0; 0; 0; dt/theta; alphaXY*dt/theta; alphaXZ*dt/theta; 0; 0; 0; 0; 0; 0; 0; 0; 0];    % for ux
B(:,2) = [0; 0; 0; 0; 0; 0; alphaXY*dt/theta; dt/theta; alphaYZ*dt/theta; 0; 0; 0; 0; 0; 0; 0; 0; 0];    % for uy
B(:,3) = [0; 0; 0; 0; 0; 0; alphaXZ*dt/theta; alphaYZ*dt/theta; dt/theta; 0; 0; 0; 0; 0; 0; 0; 0; 0];    % for uz


% Cost matrices
Qa = px_A*px_A'+py_A*py_A'+pz_A*pz_A';
Qb = px_B*px_B'+py_B*py_B'+pz_B*pz_B';
Qc = px_C*px_C'+py_C*py_C'+pz_C*pz_C';
Qn = px_C*px_C'+py_C*py_C'+pz_C*pz_C'+w_v*(vx*vx'+vy*vy'+vz*vz')+w_f*(fx*fx'+fy*fy'+fz*fz');


Q = zeros(szX,szX,N);
Q(:,:,1:N_at_A) = repmat(Qa/dim, [1,1,N_at_A]);
Q(:,:,N_at_A+1:N_at_B) = repmat(Qb/dim, [1,1,N_at_B-N_at_A]);
Q(:,:,N_at_B+1:N_at_C-1) = repmat(Qc/dim, [1,1,N_at_C-1-N_at_B]);
Q(:,:,N_at_C) = repmat(Qn/(dim*3), [1,1,1]);




%%%%% Augument matrices with feedback delay (see Crevecoeur et al., Automatica, 2011)  
nDelay = round(tDelay/dt);     % feedback delay is nDelay time steps

A = [A,zeros(szX,szX*nDelay)];
A   = [A;[eye(szX*nDelay),zeros(szX*nDelay,szX)]];

B   = [B;zeros(szX*nDelay,szU)];
H   = [zeros(szY,szX*nDelay),H];

C0  = [C0,zeros(szX,szX*nDelay)];
C0  = [C0;zeros(szX*nDelay,szX*(nDelay+1))];

E0  = [E0,zeros(szX,szX*nDelay)];
E0  = [E0;zeros(szX*nDelay,szX*(nDelay+1))];

D   = [zeros(szY,szX*nDelay),D];

Q   = cat(2,Q,zeros(szX,(szX*nDelay),N));
Q   = cat(1,Q,zeros(szX*nDelay,szX*(nDelay+1),N));

X1  = [X1;zeros(szX*nDelay,1)];
S1  = diag(repmat(diag(S1),[nDelay+1,1]));


% Scaling matrices for control dependent system noise
% The multiplicative noises affecting ux, uy and uz have the same amplitude but are uncorrelated (Crevecoeur et al., Automatica, 2011)
C(:,:,1) = B*[[sigma_u 0 0];[0 0 0];[0 0 0]];         
C(:,:,2) = B*[[0 0 0];[0 sigma_u 0];[0 0 0]];
C(:,:,3) = B*[[0 0 0];[0 0 0];[0 0 sigma_u]];

% Optimize
[K,L,~] = kalman_lqg_xstar(A,B,C,C0,H,D,D0,E0,Q,R,X1,S1);

% Internal model 
H_hat  = H;
A_hat  = A;
B_hat  = B;


xyzPos = cell(dim,1);

for d=1:dim 
    xyzPos{d,1} = zeros(Nreps,N);
end

for k=1:Nreps
    XSim = simulateOptimalReach(paramMultiplier(1)*A, paramMultiplier(2)*A_hat, paramMultiplier(3)*B, paramMultiplier(4)*B_hat, C, C0, paramMultiplier(5)*H, paramMultiplier(6)*H_hat, D, D0, E0, X1, S1, paramMultiplier(7)*L, paramMultiplier(8)*K, N, FBtype);
    switch dim
        case 1
            hold on, plot(XSim(1,:), 'Color', col);
        case 2
            hold on, plot(XSim(1,:), XSim(2,:), 'Color', col);
        case 3
            hold on, plot3(XSim(1,:), XSim(2,:), XSim(3,:), 'Color', col);
            axis equal
            set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.03 0.03], 'YDir', 'Reverse', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5, ...
        'XLim', [-1 2.5], 'YLim', [-1.3 1.3], 'ZLim', [-0.5 2], 'View', [-130 25], 'XTick', -1:0.5:2.5, 'YTick', -1:0.5:1, 'ZTick', -0.5:0.5:2)
            
    end

    for d=1:dim
        xyzPos{d,1}(k,:) = XSim(d,:);
    end
end

