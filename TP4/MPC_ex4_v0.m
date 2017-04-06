%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Model Predictive Control - Exercise 4
%              EPFL - Spring semester 2017 -
%
%            Huber Lukas - Zgraggen Jannik
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all; clear all;

%%
clc;

A = [0.9752, 1.4544; -0.0327, 0.9315];
B = [0.0248; 0.0327];

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

sys.x.penalty = QuadFunction(Q)
sys.u.penalty = QuadFunction(R)

LQRGarin = sys.LQRGain;
LQRPenalty = sys.LQRPenalty.weight;
LQRSet = sys.LQRSet;


fprintf('Programm terminated. \n')

