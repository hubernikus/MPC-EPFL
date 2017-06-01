%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Model Predictive Control - Exercise 5
%              EPFL - Spring semester 2017 - 
%
%            Huber Lukas - Zgraggen Jannik
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables;
addpath(genpath('../tbxmanager'))
addpath(genpath(' /opt/gurobi702/linux64/matlab/'))

yalmip('clear')
close all; clc;


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
%plot(refDist');
%legend('Outside Temprature in ?C','Solar gains in kW','internal gains in kW')
% JZ: write comments!!!!

d=refDist;

%% Controller Design (Setting-up MPC optimizer)

%% Section 1: tracking MPC
% MPC parameters
y_ref=[24 24 24]';         %initial conditions
umax= 15; umin=0; ymax=26; ymin=22; % constraints
R=eye(3);
Hu=[1 0 0;0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
hu=[umax,umax,umax,-umin,-umin,-umin]';
Hy=[1 0 0;0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
hy=[ymax ymax ymax -ymin -ymin -ymin]';
% Limit wrong way !!!!!

N=30;
[~,T]=size(refDist);
T=(T-N);

% Optimisation variables
x = sdpvar(10,N,'full');
y = sdpvar(3,N,'full');
u = sdpvar(3,N,'full');
d = sdpvar(3,N,'full');

% Solver options
opt = sdpsettings('verbose',1);
opt.solver = 'gurobi';

%%
% Define constraints and objective for MPC-controller
con = [];
obj = 0;

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)];
con = [con, Hu*u(:,1) <= hu];                   % Input constraints
% Include first time input constraints!!!

for j = 2:N-1  
    obj = obj + (y(:,j)-y_ref)'*R*(y(:,j)-y_ref);    % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    %/!\ d(j) or d(j+1)?
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy];                   % Output constraints
end 

%/!\ Why do we still need this constraint? Shouldn't that be fullfilled
%with [Hu*u(:,j) <= hu]; 

con=[con, u(:,:)>=0];
controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);
 

%% Section 2: economic MPC and soft constraints
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
s1 = sdpvar(6,N,'full');

% Exercise specific parameters
penal=10;
c=0.2; % CHF/kWh
R_econom=[c,c,c]; % We sample every 20min and hold the input constant over this 20 min...
 %/!\ Why not devided by 3???
% Define constraints and objective for MPC-controller

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 
% !!! input constraint yes/no?
con = [con, Hu*u(:,1) <= hu];   % Make u vector shorter?! % Input constraints

obj = R_econom*(u(:,1))


% !!! what happens with s1? how defined
for j = 2:N-1  
    obj = obj + R_econom*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)];     % Output constraints
end

con=[con, u(:,:)>=0];
%/!\ Why do we still need this constraint? Shouldn't that be fullfilled
%with [Hu*u(:,j) <= hu];  
%/!\ Why can't we put a soft constraint on u? Physical limitation?

controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);

%[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);

%% Section 3: economic, soft constraints, and variable cost
% Reset constraints and objective
con = [];
obj = 0;


% New decision varaibles
s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh

% Exercise specific parameters
penal=10; 

% Define constraints and objective for MPC-controller

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,1)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 
obj = [c(1)/3,c(1)/3,c(1)/3]*(u(:,1))

for j = 2:N-1  

    obj = obj + [c(j),c(j),c(j)]*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    %/!\ Why not devided by 3???
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)];     % Output constraints
end

con=[con, u(:,:)>=0];


controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 2);

%% Section 4 : Night setbacks
% Reset constraints and objective
con = [];
obj = 0;


% New decision varaibles
s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=10; 

% Define constraints and objective for MPC-controller

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,1)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 
obj = [c(1)/3,c(1)/3,c(1)/3]*(u(:,1))

for j = 2:N-1  

    obj = obj + [c(j),c(j),c(j)]*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    %/!\ Why not devided by 3???
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)+[sb(j);sb(j);sb(j);-sb(j);-sb(j);-sb(j)]];     % Output constraints
end

con=[con, u(:,:)>=0];


controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:);sb(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 3);

%% Section 5 : Battery coupled with the building
close all;
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
e = sdpvar(1, N,'full'); 
v = sdpvar(1, N,'full'); 
xb=sdpvar(1, N,'full');
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=10;
alpha=1;
beta=1;

% Define constraints and objective for MPC-controller

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,1)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 
obj = c(1)*e(1)

% Battery
    con = [con, v(1)==e(1)-sum(u(:,1))];
    con = [con, e(1)>=0];
    con = [con, xb(2) == xb(1)+v(1)];
    con = [con, xb(1) <= 20];
    con = [con, 0 <= xb(1)];
    con = [con, -20 <= v(1) <= 20];

for j = 2:N-1  

    obj = obj + c(j)*e(j)+  penal*s1(j)'*s1(j);   % Cost function
    %/!\ Why not devided by 3???
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)+[sb(j);sb(j);sb(j);-sb(j);-sb(j);-sb(j)]];% Output constraints
   
    % Battery
    con = [con, v(j)==e(j)-sum(u(:,j))];
    con = [con, e(j)>=0];
    con = [con, xb(j+1) == xb(j)+v(j)];
    con = [con, xb(j) <= 20];
    con = [con, 0 <= xb(j)];
   
    con = [con, -20 <= v(j) <= 20];
    
    
    
end



ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,opt,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);
[xt, yt, ut, t, et, xbt] = simBuildStorage( controller, T, @shiftPred, N);