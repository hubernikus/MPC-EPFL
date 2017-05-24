
clear;
close all;
clc;




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
hu=[umax,umax,umax,umin,umin,umin]';
Hy=[1 0 0;0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
hy=[ymax ymax ymax ymin ymin ymin]';
% Limit wrong way !!!!!

N=30;
[~,T]=size(refDist);
T=(T-N);

% Optimisation variables
x = sdpvar(10,N,'full');
y = sdpvar(3,N,'full');
u = sdpvar(3,N,'full');
d = sdpvar(3,N,'full');


%x0=x0red; %Initial condition

% Solver options
opt = sdpsettings('verbose',1);
opt.solver = 'gurobi';



% Define constraints and objective for MPC-controller
con = [];
obj = 0;

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)];
% Include first time input constraints!!!

for j = 2:N-1  
    obj = obj + (y(:,j)-y_ref)'*R*(y(:,j)-y_ref)    % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy];                   % Output constraints
end


controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);


%% Section 2: economic MPC and soft constraints
c=0.2; % CHF/kWh

s1 = sdpvar(6,N,'full');
s2 = sdpvar(6,N,'full');


% We sample every 20min and hold the input constant over this 20 min...
R_econom=diag([c/3,c/3,c/3]);

% Define constraints and objective for MPC-controller
con = [];
%obj = [(u(:,1))'*R_econom*(u(:,1))+(hu-Hu*u(:,1))'*(hu-Hu*u(:,1))+(hy-Hy*y(:,1))'*(hy-Hy*y(:,1))];

%+(hy-Hy*y(:,1))*(hy-Hy*y(:,1))
%+(hu-Hu*u(:,1))*(hu-Hu*u(:,1))

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)]; 

for j = 2:N-1  

    obj = obj + (u(:,j))'*R_econom*(u(:,j))+ s1(j)'*s1(j)+s2(j)'+s2(j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j)+ s1(j) ==hu];                   % Input constraints
    con = [con, Hy*y(:,j)+ s2(j) ==hy];  
                 % Output constraints
end

con=[con, s1(:,:)>=0];
con=[con, s2(:,:)>=0];

controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);



%% Section 3: economic, soft constraints, and variable cost

%fill in here

%% Section 4 : Night setbacks

%fill in here

%% Section 5 : Battery coupled with the building

%fill in here
