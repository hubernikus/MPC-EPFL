clc;
close all;


yalmip('clear')
clear all

%% Model data

load building.mat;
load battery.mat;

% Parameters of the Building Model
A = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C = ssM.C;
Ts = ssM.timestep;

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu;   

% Installation Test
yalmip('version')
sprintf('The Project files are successfully installed')

% Plotting the disturbance 
plotDisturbance(refDist);

%Other parameters


%% Controller Design (Setting-up MPC optimizer)

%% Section 1: tracking MPC

% Setting the parameters
N = 72; % Horizon
yRef = [24 24 24]'; % Reference
R = 1;
T = 10; % Simulation length in time-steps

% Defining the input and output constraints
ku = [15; 0];
Hu = [1; -1];

% Define optimization variables for MPC
xMPC = sdpvar(10, N,'full');
uMPC = sdpvar(3, N-1,'full');
yMPC = sdpvar(3, N,'full');
d = sdpvar(3, N,'full');
x0 = sdpvar(10,1,'full'); % Initial State

% Define the cost function for MPC
objective = 0;
for i = 1:N
    objective = objective + (yMPC(:,i)-yRef)'*R*(yMPC(:,i)-yRef);
end

% Define constraints
constraints = [];
constraints = [constraints, xMPC(:,1) == x0];
for i = 1:N-1
    % System dynamics
    constraints = [constraints, yMPC(:,i) == C*xMPC(:,i)];
    constraints = [constraints, xMPC(:,i+1) == A*xMPC(:,i) + Bu*uMPC(:,i) + Bd*d(:,i)]; 
    % Input constraints
    for j = 1:3
        constraints = [constraints, Hu*uMPC(j,i) <= ku];
    end
end
constraints = [constraints, yMPC(:,N) == C*xMPC(:,N)];

ops = sdpsettings('verbose',1);
controller = optimizer(constraints,objective,ops,[x0;d(:)],uMPC);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);


%% Section 2: economic MPC and soft constraints

% Setting the parameters
N = 72; % Horizon
p = 0.2; % Price of electricity in $/kWh
pM = [p p p]; % Price matrix
yRef = [24 24 24]'; % Reference
S = 10;
T = 300; % Simulation length in time-steps

% Defining the input and output constraints
ku = [15; 0];
Hu = [1; -1];
ky = [26; -22];
Hy = [1; -1];

% Define optimization variables for MPC
xMPC = sdpvar(10, N,'full');
uMPC = sdpvar(3, N-1,'full');
yMPC = sdpvar(3, N,'full');
epsMPC = sdpvar(3, N,'full'); % Slack variable
d = sdpvar(3, N,'full');
x0 = sdpvar(10,1,'full'); % Initial State

% Define the cost function for MPC: cost = sum(p*input)
objective = 0;
% *) minimize the cost (linked to the price)
for i = 1:N-1
    objective = objective + pM*uMPC(:,i);
end
% *) minimize the violation over the horizon
for i = 1:N
    objective = objective + epsMPC(:,i)'*S*epsMPC(:,i);
end

% Define constraints
constraints = [];
constraints = [constraints, xMPC(:,1) == x0];
for i = 1:N-1
    % System dynamics
    constraints = [constraints, yMPC(:,i) == C*xMPC(:,i)];
    constraints = [constraints, xMPC(:,i+1) == A*xMPC(:,i) + Bu*uMPC(:,i) + Bd*d(:,i)]; 
    % Input constraints
    for j = 1:3
        constraints = [constraints, Hu*uMPC(j,i) <= ku];
    end
    % Output constraints
    for j = 1:3
        constraints = [constraints, Hy*yMPC(j,i) <= ky + epsMPC(j,i)];
    end
end
% last Output constraint
for j = 1:3
    constraints = [constraints, Hy*yMPC(j,N) <= ky + epsMPC(j,N)];
end
constraints = [constraints, yMPC(:,N) == C*xMPC(:,N)];

ops = sdpsettings('solver','gurobi');
controller = optimizer(constraints,objective,ops,[x0;d(:)],uMPC);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);

%% Section 3: economic, soft constraints, and variable cost

% Setting the parameters
N = 72; % Horizon
p = 0.2; % Price of electricity in $/kWh
pM = [p p p]; % Price matrix
yRef = [24 24 24]'; % Reference
S = 10;
T = 300; % Simulation length in time-steps

% Defining the input and output constraints
ku = [15; 0];
Hu = [1; -1];
ky = [26; -22];
Hy = [1; -1];

% Define optimization variables for MPC
xMPC = sdpvar(10, N,'full');
uMPC = sdpvar(3, N-1,'full');
yMPC = sdpvar(3, N,'full');
epsMPC = sdpvar(3, N,'full'); % Slack variable
d = sdpvar(3, N,'full');
c = sdpvar(1, N,'full'); % Cost variable
x0 = sdpvar(10,1,'full'); % Initial State

% Define the cost function for MPC: cost = sum(p*input)
objective = 0;
% *) minimize the cost (linked to the price)
for i = 1:N-1
    objective = objective + [c(:,i) c(:,i) c(:,i)]*uMPC(:,i);
end
% *) minimize the violation over the horizon
for i = 1:N
    objective = objective + epsMPC(:,i)'*S*epsMPC(:,i);
end

% Define constraints
constraints = [];
constraints = [constraints, xMPC(:,1) == x0];
for i = 1:N-1
    % System dynamics
    constraints = [constraints, yMPC(:,i) == C*xMPC(:,i)];
    constraints = [constraints, xMPC(:,i+1) == A*xMPC(:,i) + Bu*uMPC(:,i) + Bd*d(:,i)]; 
    % Input constraints
    for j = 1:3
        constraints = [constraints, Hu*uMPC(j,i) <= ku];
    end
    % Output constraints
    for j = 1:3
        constraints = [constraints, Hy*yMPC(j,i) <= ky + epsMPC(j,i)];
    end
end
% last Output constraint
for j = 1:3
    constraints = [constraints, Hy*yMPC(j,N) <= ky + epsMPC(j,N)];
end
constraints = [constraints, yMPC(:,N) == C*xMPC(:,N)];

ops = sdpsettings('solver','gurobi');
controller = optimizer(constraints,objective,ops,[x0;d(:);c(:)],uMPC);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 2);

%% Section 4 : Night setbacks

% Setting the parameters
N = 72; % Horizon
p = 0.2; % Price of electricity in $/kWh
pM = [p p p]; % Price matrix
yRef = [24 24 24]'; % Reference
S = 10;
T = 300; % Simulation length in time-steps

% Defining the input and output constraints
ku = [15; 0];
Hu = [1; -1];
ky = [26; -22];
Hy = [1; -1];

% Define optimization variables for MPC
xMPC = sdpvar(10, N,'full');
uMPC = sdpvar(3, N-1,'full');
yMPC = sdpvar(3, N,'full');
epsMPC = sdpvar(3, N,'full'); % Slack variable
d = sdpvar(3, N,'full');
c = sdpvar(1, N,'full'); % Cost variable
sb = sdpvar(1, N,'full'); % Setbacks variable
x0 = sdpvar(10,1,'full'); % Initial State

% Define the cost function for MPC: cost = sum(p*input)
objective = 0;
% *) minimize the cost (linked to the price)
for i = 1:N-1
    objective = objective + [c(:,i) c(:,i) c(:,i)]*uMPC(:,i);
end
% *) minimize the violation over the horizon
for i = 1:N
    objective = objective + epsMPC(:,i)'*S*epsMPC(:,i);
end

% Define constraints
constraints = [];
constraints = [constraints, xMPC(:,1) == x0];
for i = 1:N-1
    % System dynamics
    constraints = [constraints, yMPC(:,i) == C*xMPC(:,i)];
    constraints = [constraints, xMPC(:,i+1) == A*xMPC(:,i) + Bu*uMPC(:,i) + Bd*d(:,i)]; 
    % Input constraints
    for j = 1:3
        constraints = [constraints, Hu*uMPC(j,i) <= ku];
    end
    % Output constraints
    for j = 1:3
        constraints = [constraints, Hy*yMPC(j,i) <= ky + epsMPC(j,i) + [sb(:,i); sb(:,i)]];
    end
end
% last Output constraint
for j = 1:3
    constraints = [constraints, Hy*yMPC(j,N) <= ky + epsMPC(j,N) + [sb(:,N); sb(:,N)]];
end
constraints = [constraints, yMPC(:,N) == C*xMPC(:,N)];

ops = sdpsettings('solver','+gurobi');
controller = optimizer(constraints,objective,ops,[x0;d(:);c(:);sb(:)],uMPC);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 3);

%% Section 5 : Battery coupled with the building

% Setting the parameters
N = 72; % Horizon
p = 0.2; % Price of electricity in $/kWh
pM = [p p p]; % Price matrix
yRef = [24 24 24]'; % Reference
S = 10;
T = 300; % Simulation length in time-steps
Ab = ssModel.A; %Battery state matrix
Bub = ssModel.Bu; %Battery input matrix

% Defining the input and output constraints
ku = [15; 0];
Hu = [1; -1];
ky = [26; -22];
Hy = [1; -1];
kv = [20; 20];
Hv = [1; -1];
kxb = [20; 0];
Hxb = [1; -1];

% Define optimization variables for MPC
xMPC = sdpvar(10, N,'full');
uMPC = sdpvar(3, N-1,'full'); % Power Consumption of the HVAC system
yMPC = sdpvar(3, N,'full');
epsMPC = sdpvar(3, N,'full'); % Slack variable
d = sdpvar(3, N,'full');
c = sdpvar(1, N,'full'); % Cost variable
sb = sdpvar(1, N,'full'); % Setbacks variable
x0 = sdpvar(10,1,'full'); % Initial State
e = sdpvar(1, N-1,'full'); % Power Consumption from the grid
xb = sdpvar(1, N,'full'); % State of the battery
xb0 = sdpvar(1,1,'full'); % Initial Battery State
v = sdpvar(1, N-1,'full'); % Charging power

% Define the cost function for MPC: cost = sum(p*input)
objective = 0;
% *) minimize the cost (linked to the price)
for i = 1:N-1
    objective = objective + c(:,i)*e(:,i);
end
% *) minimize the violation over the horizon
for i = 1:N
    objective = objective + epsMPC(:,i)'*S*epsMPC(:,i);
end

% Define constraints
constraints = [];
constraints = [constraints, xMPC(:,1) == x0];
constraints = [constraints, xb(:,1) == xb0];
for i = 1:N-1
    % System dynamics
    constraints = [constraints, yMPC(:,i) == C*xMPC(:,i)];
    constraints = [constraints, xMPC(:,i+1) == A*xMPC(:,i) + Bu*uMPC(:,i) + Bd*d(:,i)]; 
    % Input constraints
    for j = 1:3
        constraints = [constraints, Hu*uMPC(j,i) <= ku];
    end
    % Output constraints
    for j = 1:3
        constraints = [constraints, Hy*yMPC(j,i) <= ky + epsMPC(j,i) + [sb(:,i); sb(:,i)]];
    end
    % Battery model
    constraints = [constraints, v(:,i) == (e(:,i)-sum(uMPC(:,i),1))];
    constraints = [constraints, xb(:,i+1) == Ab*xb(:,i) + Bub*v(:,i)];
    % Battery state constraints
    constraints = [constraints, Hxb*xb(i) <= kxb];
    % Battery input constraints
    constraints = [constraints, Hv*v(:,i) <= kv];
end
% last Output constraint
for j = 1:3
    constraints = [constraints, Hy*yMPC(j,N) <= ky + epsMPC(j,N) + [sb(:,N); sb(:,N)]];
end
constraints = [constraints, yMPC(:,N) == C*xMPC(:,N)];
% last Battery state constraint
constraints = [constraints, Hxb*xb(N) <= kxb];
% we cannot reinject power into the grid
for i = 1:(N-1)
    constraints = [constraints, e(:,i) >= 0];
end


ops = sdpsettings('solver','+gurobi');
controller = optimizer(constraints,objective,ops,[x0;xb0;d(:);c(:);sb(:)],[uMPC;v;e]);
[xt, yt, ut, t, et, xbt] = simBuildStorage(controller, T, @shiftPred, N);

