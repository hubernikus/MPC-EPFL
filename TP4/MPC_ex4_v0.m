%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Model Predictive Control - Exercise 4
%              EPFL - Spring semester 2017 -
%
%            Huber Lukas - Zgraggen Yannik
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all; clear all;

%%
A = [0.9752, 1.4544; -0.0327, 0.9315];
b = [0.0248; 0.0327];

f_nat = 0.15; % [r/s] Natural frquency

discRate = 1.5; % [r/s] Discretization rate applied

xi = 0.1; % Dampin Ratio

N = 10; % Horizon length




