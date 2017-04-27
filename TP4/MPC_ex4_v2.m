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


%%
clc; 

A = [0.9752, 1.4544; -0.0327, 0.9315];   B = [0.0248; 0.0327];

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

LQRGarin = sys.LQRGain
Q_f =  sys.LQRPenalty.weight
F_f = sys.LQRSet.A
f_f = sys.LQRSet.b

% Ex2 Matlab optimization
%z = [x;u];
%H = blkdiag(kron(eye(N???1),Q), Qf, kron(eye(N),R))

% Optimization
H=blkdiag(kron(eye(N-1),Q),Q_f,kron(eye(N),R));
%H= [Q,zeros(dimX,dimU);zeros(dimU,dimX) R];
h = zeros(N*(dimX+dimU),1);

% Define Matrizes for comparison restriction
F  = [1 0; -1 0; 0 1; 0 -1]; f = [5 5 0.2 0.2]';
M = [1;-1]; m = [1.75 1.75]';

G = blkdiag(kron(eye(N-1), F), F_f,kron(eye(N), M));
g = [kron(ones(N-1,1),f);f_f; kron(ones(N,1),m)];
x0 = [3 0]'; % Starting 

% Create Equality matrizes Aeq and beq
T = [eye(N*dimX) + kron(diag(ones(1,N-1),-1),-A),  kron(diag(ones(1,N)),B)];
t = [A; zeros(dimX*(N-1),dimX)]*x0;

fprintf('Start optimization \n')
options = optimoptions(@quadprog,'ConstraintTolerance',1e0)
[zopt, fval, flag] = quadprog(H, h, G, g, T, t,[],[],[],options)

%'ConstraintTolerance',1e-8)
%                                                ConstraintTolerance            


fprintf('Programm terminated. \n')

