%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Model Predictive Control - Exercise 4
%              EPFL - Spring semester 2017 -
%
%            Huber Lukas - Zgraggen Jannik
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all; clear all;

addpath(genpath('../tbxmanager'))


%% Exercise 1
clc; 

A = [0.9752, 1.4544; -0.0327, 0.9315];   
B = [0.0248; 0.0327];
x0=[3;0];

dimX = size(A,1);
dimU = size(B,2);

Q = [10 0; 0 10];
R = [1];

f_nat = 0.15; % [r/s] Natural frquency
discRate = 1.5; % [r/s] Discretization rate applied

xi = 0.1; % Dampin Ratio

N = 10; % Horizon length


sys = LTISystem('A',A,'B',B);

% Define limits
sys.x.min = [-5, -0.2]';
sys.x.max = [5, 0.2];

sys.u.min = -1.75;
sys.u.max = 1.75;

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

LQRGarin = sys.LQRGain;
LQRPenalty = sys.LQRPenalty.weight;
LQRSet = sys.LQRSet;
Ff=LQRSet.A;
ff=LQRSet.b;
Qf=LQRPenalty;  


% Ex2 Matlab optimization
%z = [x;u];
%H = blkdiag(kron(eye(N???1),Q), Qf, kron(eye(N),R))

% Optimization
H=blkdiag(kron(eye(N-1),Q),Qf,kron(eye(N),R));
%H= [Q,zeros(dimX,dimU);zeros(dimU,dimX) R];
h = zeros(N*(dimX+dimU),1);

% Define Matrizes for comparison restriction
g = [kron(ones(N-1,1),[5 5 0.2 0.2]');ff; kron(ones(N,1),[1.75 1.75]')];
G = blkdiag(kron(eye(N-1),[1 0; -1 0; 0 1; 0 -1]),Ff, ...
            kron(eye(N),[1;-1]));

% Create Equality matrizes Aeq and beq
T = [eye(N*dimX) + kron(diag(ones(1,N-1),-1),-A),  kron(diag(ones(1,N)),B)];
t = [A; zeros(dimX*(N-1),dimX)];

t=t*x0;
options=optimoptions('quadprog','ConstraintTolerance',1e-2);

[zopt, fval, flag] = quadprog(H, h, G, g, T, t,[],[],x0,options);



% Get and plot optimal trajectory data

t_horizon=[1:N];
x1=zopt(1:2:2*(N));
x2=zopt(2:2:2*(N));
u=zopt(2*N+1:length(zopt));

h_bound=[sys.x.max;sys.x.max];
H_bound=[1 0;0 1;-1 0;0 -1 ];
P=Polyhedron([H_bound],[h_bound]);

% Check if state always lies inside boundaries
figure
P.plot
hold on
scatter(x1,x2,'b','*');

% Check if input always lies inside boundaries
figure
plot(t_horizon,u,'r');
refline(0,1.75)
refline(0,-1.75)

% Question 3: 

%/!\ How is the sytem supposed to change with u and x?

%% Exercise 2
% Parameter Definition
F = [ 1  0;             %State constraint
     -1  0;
      0  1;
      0 -1];
f = [5; 5; 0.2; 0.2];   %State constraint
m = [1.75; 1.75];       %Input constraint
M = [1; -1];            %Input constraint

% Define optimization variables
x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');

% Define constraints and objective
con = [];
obj = 0;
con = [con, x(:,1) == x0];
for i = 1:N-1      
    con = [con, x(:,i+1) == A*x(:,i) + B*u(:,i)]; % System dynamics
    con = [con, F*x(:,i) <= f];                   % State constraints
    con = [con, M*u(:,i) <= m];                   % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

%con = [con, Ff*x(:,N) <= ff];       % Terminal constraint
obj = obj + x(:,N)'*Qf*x(:,N);      % Terminal weight

% Compile the matrices
opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver
ctrl = optimizer(con, obj, opt, x(:,1), u(:,1));
% Can now compute the optimal control input using
[uopt,isfeasible] = ctrl{x0}
% isfeasible == 1 if the problem was solved successfully





fprintf('Programm terminated. \n')