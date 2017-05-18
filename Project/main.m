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


%Other parameters

%Fill in here

%% Controller Design (Setting-up MPC optimizer)

%% Section 1: tracking MPC

%fill in here

%% Section 2: economic MPC and soft constraints

%fill in here

%% Section 3: economic, soft constraints, and variable cost

%fill in here

%% Section 4 : Night setbacks

%fill in here

%% Section 5 : Battery coupled with the building

%fill in here