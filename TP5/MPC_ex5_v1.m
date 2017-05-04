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

%% System initialization
close all; clc;

A = [0.7115, -0.4345; 0.4345, 0.8853];
B = [0.2173; 0.0573];

C = [0, 1];
%D = d

% Augmented system
B_d = zeros(2,1);
C_d = [1];

A_augm = [A,B_d;zeros(1,2),1];
B_augm = [B;0];
C_augm = [C, 1];

L = (place(A_augm',-C_augm',[0.5,0.6,0.7]))'; % Why positiv?!


% Initial Estimation
x0_est = [3;0];
d_est = [0];

%Initial Conditions - Real system
x0_r = [2;1];
d_r = 0.2;



%% Exercise 1 - Observer Design
close all;

% Error
deltaX = [x0_r-x0_est];
deltaD = [d_r-d_est];

obsError = [deltaX; deltaD];

MAXITER = 40;
for i = 2:MAXITER
    obsError(:,i) = (A_augm + L *C_augm)*obsError(:,i-1);
end

figure('Position',[0 0 1000 600]); grid on;
plot(obsError(1,:)); hold on; grid on; 
plot(obsError(2,:));  
plot(obsError(3,:)); 
legend('Error x_1','Error x_2', 'Error d')

%% Exercise 2 - Steady-State Target Problem

x = sdpvar(2,1,'full');
u_r = sdpvar(1,'full');
R_s=1;
H=[1 0;0 -1];
h=[3,3]

% Define constraints and objective

%con = [con, x(:,1) == x0];
for i = 1:N-1
    con = [];
    obj = u_r*R_s*u_r;
    con = [con, [I-A,-B;C,0]*[x;u_r] == [0;0;r]]; % System dynamics
    con = [con, H*x <= h];                   % Input constraints
    
    opt = sdpsettings;
    opt.solver = 'quadprog';
    opt.quadprog.TolCon = 1e-16;

    %opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver
    ctrl = optimizer(con, obj, opt, x, u_r);

 
end
% Compile the matrices
opt = sdpsettings;
opt.solver = 'quadprog';
opt.quadprog.TolCon = 1e-16;

 %opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver
ctrl = optimizer(con, obj, opt, x(:,1), u_r(:,1));



%% Exercise 3 - Controller Design
close all;
%%%% TODO! P etc.


% Constraints
m = [3; 3];             %Input constraint
M = [1; -1];            %Input constraint


% Weights Controler that is able to track constant output reference
Q = eye(size(A,1));
R = 1;

N =5; % Horizon length

% LQR gain
sys = LTISystem('A',A,'B',B);

sys.u.min = -m(1);
sys.u.max = m(2);

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

LQRGarin = sys.LQRGain;
LQRPenalty = sys.LQRPenalty.weight;
LQRSet = sys.LQRSet;
Ff=LQRSet.A;
ff=LQRSet.b;
Qf=LQRPenalty;  

% for i=2:MAXITER
%    %u_r(i)= -K*y_r(i-1);
%    x0_r(:,i)= A*x0_r(:,i-1)+B(i)*u_r(i);
%    y_r(i)=C*x0_r(:,i)+d_r;
%    
% end

% Define optimization variables
x = sdpvar(2,N,'full');
u_r = sdpvar(1,N,'full');

% Define constraints and objective
con = [];
obj = 0;
%con = [con, x(:,1) == x0];
for i = 1:N-1      
    con = [con, x(:,i+1) == A*x(:,i) + B*u_r(:,i)]; % System dynamics
    con = [con, M*u_r(:,i) <= m];                   % Input constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u_r(:,i)'*R*u_r(:,i); % Cost function
end
obj = obj + x(:,N)'*Qf*x(:,N);      % Terminal weight

% Compile the matrices
opt = sdpsettings;
opt.solver = 'quadprog';
opt.quadprog.TolCon = 1e-16;

 %opt = sdpsettings('solver','sedumi','verbose',0); % choosing the solver
ctrl = optimizer(con, obj, opt, x(:,1), u_r(:,1));


% Initial conditions
xi = x0_est; % try to controll estimation

xd_est = [x0_est; d_est]
y_est = [C*xd_est(1:2,1)+d_est]; 

% Real conditions
y_r = [0];
u_r = [0];

x_r = [x0_r];
y_r = [C*xi+d_r];


t = [0];

MAXITER = 100; tolX = 1e-8;
% Can now compute the optimal control input using
for i = 2:MAXITER
    [u_opt,infeasible] = ctrl{xi};
    
    if(infeasible); fprintf('Problem infeasible at i=%d \n',i); break; end;
    
    u_r(i) = u_opt;
    t(i) = i;
    
    % Real sytem
    x_r(:,i) = A*x_r(:,i-1) + B*u_opt;
    y_r(i) = C*x_r(:,1)+d_r;
        
    % Estimated system
%    xdA= A_augm*xd_est(:,i-1)
%    xdB = B_augm*u_opt
%    xdL = L*(C*xd_est(1:2,i-1) + C_d*u_opt - y_r(i-1))
    xd_est(:,i) = [A_augm*xd_est(:,i-1)+B_augm*u_opt ...
                  + L*(C*xd_est(1:2,i-1) + C_d*u_opt - y_r(i-1))];

    %y_est(i) = C*xd_est(1:2,i)+d_est
    
    xi = xd_est(1:2,i);
                
        
    if(norm(xd_est(1:2,i)-xd_est(1:2,i-1)) < tolX);
        fprintf('System converged at after %d steps. \n',i);
        break
    end
end

%x_val = [x_val, xi]; % save current values 
%t_val = [t_val, i];

%
close all;

figure('Position',[0 0 1000 600]);
plot(xd_est(1,:),xd_est(2,:),'b-*');
grid on; hold on;
plot(x_r(1,:),x_r(2,:),'-*');
legend('Estimation','Real System')
xlabel('x_1'), ylabel('x_2')        


figure('Position',[0 0 1000 600]); grid on;
plot(u_r); hold on; grid on; 
plot(y_r);
legend('u(t)','y(t)')
        


%%
close all;
%for i=2:MAXITER
%    %u(i) = 
%     x0_est(:,i)= A_augm*x0_est(:,i-1)+L*((C*x0_est([1 2],i-1))+C_d*x0_est(3,i-1)-y(i-1));
%    y(i) = C*x0_est([1,2],i)+C_d*x0_est(3,i);
%    
% end

% figure('Position',[0,0,600,800])
% hold on
% plot(x0_r(1,:),x0_r(2,:),'-*');
% plot(x0_est(1,:),x0_est(2,:),'-o');
% xlabel('x_1'),ylabel('x_2')
% legend('System','Estimation')
% print('EstimationX1X2','-dpng')
% 
% figure('Position',[0,0,600,800])
% plot([1:length(y)],y,[1:length(y)],y_r)
% xlabel('Step n'), ylabel('y')
% legend('y','y_r')
% % We did not get rid of the disturbance



%% 
fprintf('Programm terminated. \n')