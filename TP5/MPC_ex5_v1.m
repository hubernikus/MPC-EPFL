%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Model Predictive Control - Exercise 5
%              EPFL - Spring semester 2017 - 
%
%            Huber Lukas - Zgraggen Jannik
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; close all; clear all;

addpath(genpath('../tbxmanager'))

%% System initialization
close all;
clear all;
A = [0.7115, -0.4345; 0.4345, 0.8853];
B = [0.2173; 0.0573];

C = [0, 1];
%D = d
B_d=zeros(2,1);
C_d=[1];
A_er=[A,B_d;zeros(1,2),1];
C_er=[C, 1];
L=(place(A_er',-C_er',[0.1,0.2,0.3]))';
% System dynamics with estimated ..
x=[3;0;0];
y=[0];
for i=2:50
   
   x(:,i)= A_er*x(:,i-1)+L*((C_er*x(:,i-1))-y(i-1));
   y(i)=C*x([1,2],i)+x(3,i);
    
end

%True system dynamics
x_r=[3;0];
d_r=12000000;
y_r=[0];

for i=2:50
   x_r(:,i)= A*x_r(:,i-1);
   y_r(i)=C*x_r(:,i)+d_r;
    
end
figure
hold on
scatter(x_r(1,:),x_r(2,:),'*');
scatter(x(1,:),x(2,:));


%% Exercise 1 - Observer Design
x0 = [1;2];
x0_hat = [3;0]; % Estimate of x0

d = 0.2;
d_hat = 0; % Estimate of d

u = 0;

N = 10; % Number of prediction steps

% Parameter Definition
m = [3; 3];       %Input constraint
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
%    con = [con, F*x(:,i) <= f];                   % State constraints
    con = [con, M*u(:,i) <= m];                   % Input constraints
%    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

% con = [con, Ff*x(:,N) <= ff];       % Terminal constraint
% obj = obj + x(:,N)'*Qf*x(:,N);      % Terminal weight

% Compile the matrices
opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver
obs = optimizer(con, obj, opt, x(:,1), u(:,1));
% Can now compute the optimal control input using
[uopt,isfeasible] = ctrl{x0}
% isfeasible == 1 if the problem was solved successfully



%% Exercise 2 - 


