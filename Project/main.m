<<<<<<< HEAD
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
=======

clear;
close all;
clc;




clc;
close all;
>>>>>>> 883f25dc072c0a4d4a9919f4c9b6b57bfc01471d

%% System initialization

A = [0.7115, -0.4345; 0.4345, 0.8853];
B = [0.2173; 0.0573];

C = [0, 1];

% Augmented system
B_d = zeros(2,1);
C_d = [1];

A_augm = [A,B_d;zeros(1,2),1];
B_augm = [B;0];
C_augm = [C, 1];


% Make sure that eigenvalues of (A+LC) are in unit circle
L = (place(A_augm',-C_augm',[0.5,0.6,0.7]))'; 

% Initial Estimation
x0_est = [3;0];
d0_est = [0];

%Initial Conditions - Real system
x0_r = [1;2];
d_r = 0.2;

%% Exercise 1 - Observer Design

% Error between real system and estimation
deltaX = [x0_r-x0_est];
deltaD = [d_r-d0_est];

obsError = [deltaX; deltaD];

<<<<<<< HEAD

% Rund the integral disturbance dynamics
=======
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
>>>>>>> 883f25dc072c0a4d4a9919f4c9b6b57bfc01471d

MAXITER = 50; minTol = 1e-2;

<<<<<<< HEAD
for i = 2:MAXITER
    obsError(:,i) = (A_augm + L *C_augm)*obsError(:,i-1);
    if(norm(obsError(i)) < minTol)
        fprintf('Problem converged after iteration %d \n',i);
        break;
    end
end

% Plot the results

figure('Position',[0 0 1000 600]); grid on;
plot(sqrt(sum(obsError(1:2,:).^2,1))); hold on;
plot(obsError(3,:)); grid on;
legend('Error x', 'Error d')
xlabel('Step [s]'); ylabel('Error')
title('Estimator error');


% Estimation converges very nicely towards the real value. The error in x
% and d converges below the minimum Tolerance in less than 100 iterations.

% Initialize vectores
xVal_est = [x0_est]; xVal_r = [x0_r];
dVal_est = [d0_est]; %dVal_r = [d_r];
    

% Define loop
MAXITER= MAXITER; minTol = 1e-2;

for i = 2:MAXITER
    u = 0; % No control or input
    y = C*xVal_r(:,i-1)+d_r;
    
    xVal_r(:,i) = A*xVal_r(:,i-1) + B*u;
    
    dh_hat2 = A_augm*[xVal_est(:,i-1); dVal_est(i-1)]+B_augm*u ...
            +  L*(C*xVal_est(:,i-1)+C_d*dVal_est(i-1)-y);
    xVal_est(:,i) = dh_hat2(1:2);
    dVal_est(i) = dh_hat2(3);
    
    if(and(norm(xVal_est(i)-xVal_r(i)) < minTol, ...
            norm(xVal_r(:,i)-xVal_r(:,i-1)) < minTol))
        fprintf('Problem converged after iteration %d \n',i);
        break;
    end
end
    

figure('Position',[0 0 1000 600]); hold on; grid on;
plot(xVal_est(1,:), xVal_est(2,:),'b-*')
plot(xVal_r(1,:),xVal_r(2,:), 'r-*')
xlabel('x_1'); ylabel('x_2')
legend('Estimation','Real System')
title('Comparison between estimated and real state');

figure('Position',[0 0 1000 600]); hold on; grid on;
plot(dVal_est,'b'); plot([0, length(dVal_est)], [d_r, d_r], 'r')
xlabel('Step'); ylabel('Error')
legend('Real disturbance', 'Estimated disturbance')
title('Comparison of real and estimated disturbance');



%% Exercise 2 & 3 - Controller Design


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


% Weight of final cost
P = dlyap(A,Q);

% Solver settings
opt = sdpsettings;
opt.solver = 'quadprog';
opt.quadprog.TolCon = 1e-16;


% Initial conditions
xi = x0_est; % try to controll estimation
xd_est = [x0_est; d0_est]
y_est = [C*xd_est(1:2,1)+d0_est]; 
r_val = [0.5, 1];

for r = r_val
    % Real conditions
    y_r = [0];
    u_r = [0];

    x_r = [x0_r];
    y_r = [C*xi+d_r];

    t = [0];

    tolX = 1e-8;
    % Can now compute the optimal control input using
    for i = 2:MAXITER
        
        % Exercise 2 - Optimize u^2 (Target tracking)
        x_s = sdpvar(2,1,'full');
        u_s = sdpvar(1,1,'full');
        obj_ss = u_s*R*u_s;
        con_ss =[I-A,-B;C,0]*[x_s;u_s] == [0;0;r-C_d*xd_est(3,i-1)]; % System dynamics
        con_ss = [con_ss, H*u_s <= h];                       % Input constraint
        solvesdp(con_ss, obj_ss, opt);
        x_s=double(x_s);
        u_s=double(u_s);

        % Define constraints and objective for MPC-controller
        con = [];
        obj = 0;
        %con = [con, x(:,1) == x0];
        for j = 1:N-1  
            obj = obj + (x(:,j)-x_s)'*Q*(x(:,j)-x_s) + (u(:,j)-u_s)'*R*(u(:,j)-u_s); % Cost function
            con = [con, x(:,j+1) == A*x(:,j) + B*u(:,j)]; % System dynamics
            con = [con, H*u(:,j) <= h];                   % Input constraints
        end
        obj = obj + (x(:,N)-x_s)'*P*(x(:,N)-x_s);      % Terminal weight
        ctrl = optimizer(con, obj, opt, x(:,1), u(:,1));
        [u_opt,infeasible] = ctrl{xi};

        if(infeasible); fprintf('Problem infeasible at i=%d \n',i); break; end;

        u_r(i) = u_opt;
        t(i) = i;

        % Real sytem
        x_r(:,i) = A*x_r(:,i-1) + B*u_opt;
        y_r(i) = C*x_r(:,i)+d_r;

        % Estimated system
        xd_est(:,i) = [A_augm*xd_est(:,i-1)+B_augm*u_opt ...
                      + L*(C*xd_est(1:2,i-1) + C_d*xd_est(3,i-1) - y_r(i-1))];

        xi = xd_est(1:2,i);

        if(norm(abs(y_r(i)-r)) < tolX);
            fprintf('System converged at after %d steps. \n',i);
            break
        end

    end
    
    % Plot results
    figure('Position',[0 0 1000 600]);
    plot(xd_est(1,:),xd_est(2,:),'b-*');
    grid on; hold on;
    plot(x_r(1,:),x_r(2,:),'r-*');
    legend('Estimation','Real System')
    xlabel('x_1'), ylabel('x_2')
    title(['Controller performance for r=',num2str(r)]);


    figure('Position',[0 0 1000 600]); grid on;
    plot(0:length(u_r)-1, u_r,'b'); hold on; grid on; 
    plot(0:length(u_r)-1,y_r,'r');
    plot([0,length(u_r)-1],[r,r],'r--')
    plot([0,length(u_r)-1],[-3,-3],'k--')
    plot([0,length(u_r)-1],[3,3],'k--')
    xlim([0,length(u_r)-1])
    ylim([-3.1,3.1])
    title(['Controller performance for r=',num2str(r)]);
    
    legend('u(t)','y(t)','reference','input constraints')
    xlabel('Steps')
       
    

end 
    

%% 
fprintf('Programm terminated. \n')
=======
% Solver options
opt = sdpsettings('verbose',1);
opt.solver = 'gurobi';

%%
% Define constraints and objective for MPC-controller
con = [];
obj = 0;

con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)];
% Include first time input constraints!!!

for j = 2:N-1  
    obj = obj + (y(:,j)-y_ref)'*R*(y(:,j)-y_ref)    % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j+1)]; % System dynamics
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
obj = R_econom*(u(:,1))

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

>>>>>>> 883f25dc072c0a4d4a9919f4c9b6b57bfc01471d
