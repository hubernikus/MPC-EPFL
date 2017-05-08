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

% Make sure that eigenvalues of (A+LC) are in unit circle
L = (place(A_augm',-C_augm',[0.5,0.6,0.7]))'; 
%/!\ Why positiv?!?

% Initial Estimation
x0_est = [3;0];
d_est = [0];

%Initial Conditions - Real system
x0_r = [1;2];
d_r = 0.2;

%% Exercise 1 - Observer Design

% Error between real system and estimation
deltaX = [x0_r-x0_est];
deltaD = [d_r-d_est];

obsError = [deltaX; deltaD];

% Rund the integral disturbance dynamics
MAXITER = 40;
for i = 2:MAXITER
    obsError(:,i) = (A_augm + L *C_augm)*obsError(:,i-1);
end

% Plot the results
figure('Position',[0 0 1000 600]); grid on;
plot(obsError(1,:)); hold on; grid on; 
plot(obsError(2,:));  
plot(obsError(3,:)); 
legend('Error x_1','Error x_2', 'Error d')


%% Exercise 2 & 3 - Controller Design
close all;

N =5; % Horizon length

% Define optimization variables
x = sdpvar(2,N,'full');
u = sdpvar(1,N,'full');

% Constraints
h = [3; 3];             %Input constraint
H = [1; -1];            %Input constraint

% Stage cost
%Weights Controller that is able to track constant output reference
Q = eye(size(A,1));
R = 1;
I=eye(2);

% Initial conditions
xi = x0_est; % try to controll estimation
xd_est = [x0_est; d_est];
y_est = [C*xd_est(1:2,1)+d_est]; 
r=0.5;

% Real conditions
y_r = [0];
u_r = [0];
x_r = [x0_r];
y_r = [C*xi+d_r];

% Weight of final cost
P = dlyap(A,Q);

% Solver settings
    opt = sdpsettings;
    opt.solver = 'quadprog';
    opt.quadprog.TolCon = 1e-16;

% Time counter
t = [0];

% Calculate the optimal control
MAXITER = 50; tolX = 1e-8;
for i = 2:MAXITER
    
    % Offset-free tracking
    x_s = sdpvar(2,1,'full');                   %Definition Variables
    u_s = sdpvar(1,1,'full');
    obj_ss = u_s*R*u_s;                         %Objective Function
    con_ss =[I-A,-B;C,0]*[x_s;u_s] ...          %System dynamics
        == [0;0;r-C_d*xd_est(3,i-1)]; 
    con_ss = [con_ss, H*u_s <= h];              %Input constraint
    solvesdp(con_ss, obj_ss, opt);
    x_s=double(x_s);                            %Output
    u_s=double(u_s);
    
    % MPC Regulator
    con = [];
    obj = 0;
    %con = [con, x(:,1) == x0];
    for j = 1:N-1  
        obj = obj + (x(:,j)-x_s)'*Q*(x(:,j)-x_s) + (u(:,j)-u_s)'*R*(u(:,j)-u_s); % Cost function

        con = [con, x(:,j+1) == A*x(:,j) + B*u(:,j)]; %System dynamics
        con = [con, H*u(:,j) <= h];                   %Input constraints
    end
    obj = obj + (x(:,N)-x_s)'*P*(x(:,N)-x_s);         %Terminal weight
    ctrl = optimizer(con, obj, opt, x(:,1), u(:,1));  %Optimal poliy
    [u_opt,infeasible] = ctrl{xi};
    
    % Check for in feasability

    if(infeasible); fprintf('Problem infeasible at i=%d \n',i); break; end;
    
    u_r(i) = u_opt;
    t(i) = i;
    
    
    x_r(:,i) = A*x_r(:,i-1) + B*u_opt;          %Real sytem
    y_r(i) = C*x_r(:,i)+d_r;
        
    % Estimated system
%    xdA= A_augm*xd_est(:,i-1)
%    xdB = B_augm*u_opt
%    xdL = L*(C*xd_est(1:2,i-1) + C_d*u_opt - y_r(i-1))
    
    % Get the estimation
    xd_est(:,i) = [A_augm*xd_est(:,i-1)+B_augm*u_opt ...
                  + L*(C*xd_est(1:2,i-1) + C_d*xd_est(3,i-1) - y_r(i-1))];

    %y_est(i) = C*xd_est(1:2,i)+d_est
    
    xi = xd_est(1:2,i);
                
    % Check for convergence    
    if(norm(abs(y_r(i)-r)) < tolX);
        fprintf('System converged at after %d steps. \n',i);
        break
    end
end

%x_val = [x_val, xi]; % save current values 
%t_val = [t_val, i];

%%
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
fprintf('Programm terminated. \n')