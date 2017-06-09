clear variables;
addpath(genpath('../tbxmanager'))

clc;

yalmip('clear')

%% Model data
close all;

load building.mat;
load battery.mat;

% Parameters of the Building Model
A  = ssM.A;
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
timeRange = 1/3*(0:size(refDist,2)-1);

figure('Position',[0 0 2000 ]
subplot(3,1,1)
plot(timeRange,refDist(1,:)); hold on; grid on;
xlim([timeRange(1), timeRange(end)]); 
ylabel('Degree C')
legend('Outside Temprature')
subplot(3,1,2)
plot(timeRange,refDist(2,:)); hold on; grid on;
legend('Solar gains')
xlim([timeRange(1), timeRange(end)]); 
ylabel('kW')
subplot(3,1,3)
plot(timeRange,refDist(3,:)); hold on; grid on;
xlabel('Time [h]')
xlim([timeRange(1), timeRange(end)]); 
ylabel('kW')
legend('internal gains')

print('fig/errorPlot','-dsvg')
set(gcf, 'PaperUnits', 'centimeters');
x_width=10;y_width=20
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,'fig/MPC_error.svg')

% JZ: write comments!!!!

d=refDist;

%% Controller Design (Setting-up MPC optimizer)

%% Section 1: tracking MPC
close all;

% MPC parameters
T_s = 20; %[min] - Sampling Time1

y_ref=[24 24 24]';         % reference conditions

umin=0;  %input constrraints
umax= 15; 

ymin=22; %  output constraints
ymax=26; 

R=eye(3);

Hu = [eye(3);eye(3)];
hu = [-1*umin*ones(3,1); umax*ones(3,1)];
%Hu_max = eye(3);
%hu_max = umin*ones(3,1);

Hy=[-1*eye(3);eye(3)];
hy= [-ymin*ones(3,1);ymax*ones(3,1)];
%Hy_max=eye(3);
%hy_max=ymax*ones(3,1);

N=50; % Number of prediction steps

[jlj,T]=size(refDist);

T=T-N;
T = 20
% Optimisation variables
x = sdpvar(10,N,'full'); % Ten states with any physical representation
y = sdpvar(3,N,'full');
u = sdpvar(3,N,'full');
d = sdpvar(3,N,'full');

%x0=x0red; %Initial condition


% Define constraints and objective for MPC-controller
con = [];
obj = 0;

%con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
%con = [con, y(:,1) == C*x(:,1)]; 
% why ?! ^^
for j = 1:N-1  
    obj = obj + (y(:,j)-y_ref)'*R*(y(:,j)-y_ref);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    %con = [con, Hu_max*u(:,j) <= hu_max];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy];                   % Output constraints
    %con = [con, Hy_max*y(:,j) <= hy_max];                   % Output constraints
end

% Solver options
opt = sdpsettings;
opt.solver = 'quadprog';
opt.quadprog.TolCon = 1e-16;
%opt = sdpsettings('verbose',1);
controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1);


%% Test

close all;

A = [0.9752, 1.4544; -0.0327, 0.9315];   
B = [0.0248; 0.0327];
x0=[3;0];

dimX = size(A,1);
dimU = size(B,2);

Q = 10*[1 0; 0 1];
R = [1];

f_nat = 0.15; % [r/s] Natural frquency
discRate = 1.5; % [r/s] Discretization rate applied

dx_i = 0.1; % Dampin Ratio

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

for i = 1:N-1     
     con = [con, x(:,i+1) == A*x(:,i) + B*u(:,i)]; % System dynamics
     con = [con, F*x(:,i) <= f];                   % State constraints
     con = [con, M*u(:,i) <= m];                   % Input constraints
     obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i); % Cost function
end

%con = [con, Ff*x(:,N) <= ff];       % Terminal constraint
%obj = obj + x(:,N)'*Qf*x(:,N);      % Terminal weight

% Compile the matrices
%opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver

opt = sdpsettings;
opt.solver = 'gurobi';
opt.quadprog.TolCon = 1e-16;
ctrl = optimizer(con, obj, opt, x(:,1), u(:,1));

% Can now compute the optimal control input using
xi = x0;
x_yalm = [];
u_yalm = [];
t_yalm = [];

infeasible = 0; i = 1;

maxIter = 100;
while(not(infeasible==1) && sum(Ff*xi>ff) > 0)
    [uOpt,infeasible] = ctrl{xi};
    
    yalmiperror(infeasible)
    x_yalm = [x_yalm, xi]; % save current values 
    u_yalm = [u_yalm, uOpt];
    t_yalm = [t_yalm, i];
    
    xi = A*xi  + B*uOpt; % Update step
    
    if(i > maxIter); fprintf('Maximum Iteration reached i=%d \n',i); break; end;
    
    i = i + 1;
    
end


%% Section 2: economic MPC and soft constraints







%% Section 3: economic, soft constraints, and variable cost

%fill in here

%% Section 4 : Night setbacks

%fill in here

%% Section 5 : Battery coupled with the building

%fill in here