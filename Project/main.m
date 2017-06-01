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
N=30;
[~,T]=size(refDist);
T=(T-N);

% Cost calculation
% Variable Cost
Thours=T/3;
Tdays=Thours/24;
refCost =0.2*ones(1,length(refDist));
for i=0:floor(Tdays)  
refCost(i*24*3+(30:30+18))=0.04;
end
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
%controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
%[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);
 

%% Section 2: economic MPC and soft constraints
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
s1 = sdpvar(6,N,'full');

% Exercise specific parameters
penal=10;
c=0.2; % CHF/kWh
R_econom=[c/3,c/3,c/3]; % We sample every 20min and hold the input constant over this 20 min...
 %/!\ Why not devided by 3???
% Define constraints and objective for MPC-controller

% Initial constraints and obj // Without constraint on y
obj = R_econom*(u(:,1))
con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 

% !!! what happens with s1? how defined
for j = 2:N-1  
    obj = obj + R_econom*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)];     % Output constraints
end
% Final constraints and obj 
con = [con, Hu*u(:,N) <= hu];
con = [con, Hy*y(:,N) <= hy + s1(N)]


controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);


Total_sc=sum(ut(:))*c/3;
%% Section 3: economic, soft constraints, and variable cost
% Reset constraints and objective
con = [];
obj = 0;


% New decision varaibles
s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh

% Exercise specific parameters
penal=1; 

% Initial constraints and obj // Without constraint on y
 
for j = 1:N-1  
    obj = obj + [c(j)/3,c(j)/3,c(j)/3]*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)];     % Output constraints
end

% Final constr
obj = obj + [c(N)/3,c(N)/3,c(N)/3]*(u(:,N))+  penal*s1(N)'*s1(N);   % Cost function
con = [con, y(:,N) == C*x(:,N)];
con = [con, Hu*u(:,N) <= hu];                   % Input constraints
con = [con, Hy*y(:,N) <= hy + s1(N)];     % Output constraints

ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 2);

Total_vc=refCost(1:T)/3*ut(1,:)'+refCost(1:T)/3*ut(2,:)'+refCost(1:T)/3*ut(3,:)';

%% Section 4 : Night setbacks
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=1; 
% Define constraints and objective for MPC-controller

for j = 1:N-1  
    obj = obj + [c(j)/3,c(j)/3,c(j)/3]*(u(:,j))+  penal*s1(j)'*s1(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];     % Output constraints
end

% Final constraints
    obj = obj + [c(N)/3,c(N)/3,c(N)/3]*(u(:,N))+  penal*s1(N)'*s1(N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + s1(N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];     % Output constraints

ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:);sb(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 3);


Total_sb=refCost(1:T)/3*ut(1,:)'+refCost(1:T)/3*ut(2,:)'+refCost(1:T)/3*ut(3,:)';

%% Section 5 : Battery coupled with the building

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

for j = 1:N-1  
    % System
    obj = obj + c(j)*e(j)+  penal*s1(j)'*s1(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + s1(j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
   
    % Battery
    con = [con, v(j)==e(j)-sum(u(:,j))];
    con = [con, e(j)>=0];
    con = [con, xb(j+1) == xb(j)+v(j)];
    con = [con, xb(j) <= 20];
    con = [con, 0 <= xb(j)];
    con = [con, -20 <= v(j) <= 20];  
end

% Final constraints
    % System
    obj = obj + c(N)*e(N)+  penal*s1(N)'*s1(N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + s1(N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
    % Battery
    con = [con, v(N)==e(N)-sum(u(:,N))];
    con = [con, e(N)>=0];
    con = [con, xb(N) <= 20];
    con = [con, 0 <= xb(N)];
    con = [con, -20 <= v(N) <= 20];  


ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,opt,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);

[xt, yt, ut, t, et, xbt] = simBuildStorage( controller, T, @shiftPred, N);

Total_bat=refCost(1:T)/3*et(1,:);
