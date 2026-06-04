dt = 0.01;                  % time step (s)
tDelay = 0.04;              % sensory delay (s)
m = 1;                      % mass 
b = 0.1;                    % viscous damping coefficient.
theta = 0.03;               % time constant (s)
p_zero = [0 0 0];           % [x y z]
p_trgA = [1.6 0.25 1.3];    % [x y z]
p_trgB = [1.6 -0.5 1.1];    % [x y z]
p_trgC = [0 0 0.2];         % [x y z]

T_end = 0.53;               % movement duration (s)
T_A = 0.2;                  % reach duration (s)
T_B = 0.08;                 % grasp duration (s)
T_C = 0.25;                 % return duration (s)

       
% coupling between actions in different dimensions (see Crevecoeur et al., Automatica, 2011)   
alphaXY = 0;      
alphaXZ = 0;
alphaYZ = 0;

sigma_u = 0.1;              % control noise scale factor 
sigma_s = 0.3;              % sensory noise scale factor

dim = 3;                    % movement dimension

Wt = [0.01; 0.1; 0];        % vector of state-dependent sensory noise terms: [position; velocity; force]

w_v = 1e-1;                 % scaling weight for the velocity cost 
w_f = 1e-2;                 % scaling weight for the force cost
r = 5e-6;                   % weight of the energy term

S1scale = 0.002;
E0scale = 0;
D0scale = 0.5;
C0scale = 0.01;

FBtype = 'Crevecoeur2011';

Nreps = 100;                % number of simulated trajectories