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
LQRPenalty = sys.LQRPenalty.weight
LQRSet = sys.LQRSet


% Ex2 Matlab optimization
% z = [x;u];
%H = blkdiag(kron(eye(N−1),Q), Qf, kron(eye(N),R))

% Optimization
H= [Q,zeros(dimX,dimU);zeros(dimU,dimX) R];
h = zeros(dimX+dimU,1);
%H = blkdiag(kron(eye(N−1),Q), Qf, kron(eye(N),R))

% Define Matrizes for comparison restriction
%g = kron([-1;15 5 0.2 ]');
g = [kron(ones(N,1),[5 -5 0.2 -0.2]'); kron(ones(N,1),[1.75 -1.75]')];

G = blkdiag(kron(eye(N),[1 0; -1 0; 0 1; 0 -1]), ...
            kron(eye(N),[1;-1]))


% Create Equality matrizes Aeq and beq
T = [eye(N*dimX) + kron(diag(ones(1,N-1),-1),A),  kron(diag(ones(1,N)),B)];
t = [A; zeros(dimX*(N-1),dimX)];

%[zopt, fval, flag] = quadprog(H, h, G, g, T, t);




!!!! REMOV
fprintf('Programm terminated. \n')

