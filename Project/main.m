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
addpath(genpath(' /../../../../../../../../opt/gurobi702/linux64/matlab/'))

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

%Fill in here
timeRange = 1/3*(0:size(refDist,2)-1);

figure('Position',[0 0 1000 400])
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
print('fig/MPC_externalInput','-dpng')


d=refDist;
N=72;
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
%% Section 1: tracking MPC

% MPC parameters
y_ref=[24 24 24]';         %initial conditions
umax= 15; umin=0; ymax=26; ymin=22; % constraints
%R=eye(3);
R=diag([1 2 6]);
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

% Constraints and Objective
con = [];
obj = 0;
con = [con, x(:,2) == A*x(:,1) + Bu*u(:,1)+Bd*d(:,2)]; % System dynamics
con = [con, y(:,1) == C*x(:,1)];
con = [con, Hu*u(:,1) <= hu];                   % Input constraints

for j = 2:N-1  
    obj = obj + (y(:,j)-y_ref)'*R*(y(:,j)-y_ref);    % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy];                   % Output constraints
end 

obj = obj + (y(:,N)-y_ref)'*R*(y(:,N)-y_ref);
con = [con, y(:,N) == C*x(:,N)];
con = [con, Hu*u(:,N) <= hu]; 
con = [con, Hy*y(:,N) <= hy];  

% Simulation
controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
[xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1,'fig/firstMPC');
 
totE_MPC = sum(ut);
totC_MPC = 0.3*0.2*sum(ut);


%% Section 2: economic MPC and soft constraints
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
S = sdpvar(6,N,'full');

% Exercise specific parameters

penal=1;

c=0.2; % CHF/kWh
R_econom=[c/3,c/3,c/3]; % We sample every 20min and hold the input constant over this 20 min...

% Define constraints and objective for MPC-controller

% Constraints and Objectives
for j = 1:N-1  
    obj = obj + R_econom*(u(:,j))+  penal*S(:,j)'*S(:,j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)];     % Output constraints
    saveS(:,j) = S(:,j);
end
    obj = obj + R_econom*(u(:,N))+  penal*S(:,N)'*S(:,N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + S(:,N)];     % Output constraints
  
% % Simulation
% opt = sdpsettings('verbose',1, 'solver', '+gurobi');
% controller = optimizer(con,obj,opt,[x(:,1);d(:)],u);
% 
% [xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 1,'fig/softConstrK');
% 
% % Total Cost
% totCost_sc=sum(ut(:))*c/3;
% totEner_sc = sum(ut);


%% Section 3: economic, soft constraints, and variable cost
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
c = sdpvar(1, N,'full'); % CHF/kWh

% Exercise specific parameters
penal=1; 

% Constraints and Objectives
for j = 1:N-1  
    obj = obj + [c(j)/3,c(j)/3,c(j)/3]*(u(:,j))+  penal*S(:,j)'*S(:,j); % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)];     % Output constraints
end
obj = obj + [c(N)/3,c(N)/3,c(N)/3]*(u(:,N))+penal*S(:,N)'*S(:,N);   % Cost function
con = [con, y(:,N) == C*x(:,N)];
con = [con, Hu*u(:,N) <= hu];                   % Input constraints
con = [con, Hy*y(:,N) <= hy + S(:,N)];     % Output constraints

% % Simulation
% opt = sdpsettings('verbose',1, 'solver', '+gurobi');
% controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:)],u);
% [xt, yt, ut, t] = simBuild(controller, T, @shiftPred, N, 2,'fig/varCost');
% 
% % Total Cost
% totCost_vc=refCost(1:T)/3*ut(1,:)'+refCost(1:T)/3*ut(2,:)'+refCost(1:T)/3*ut(3,:)';
% totEner_vc = sum(ut);

%% Section 4 : Night setbacks
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles

c = sdpvar(1, N,'full'); % CHF/kWh
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=1; 
% Define constraints and objective for MPC-controller

for j = 1:N-1  
    obj = obj + [c(j)/3,c(j)/3,c(j)/3]*(u(:,j))+  penal*S(:,j)'*S(:,j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];     % Output constraints
end

% Final constraints
    obj = obj + [c(N)/3,c(N)/3,c(N)/3]*(u(:,N))+  penal*S(:,N)'*S(:,N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];     % Output constraints

% % Simulation
% 
% ops = sdpsettings('verbose',1, 'solver', '+gurobi');
% controller = optimizer(con,obj,opt,[x(:,1);d(:);c(:);sb(:)],u);
% [xt_sb, yt_sb, ut_sb, t_sb] = simBuild(controller, T, @shiftPred, N, 3,'fig/nightTime');

% et_sb=ut_sb(1,:)+ut_sb(2,:)+ut_sb(3,:);
% Total_sb=refCost(1:T)/3*ut_sb(1,:)'+refCost(1:T)/3*ut_sb(2,:)'+refCost(1:T)/3*ut_sb(3,:)';

%% Section 5 : Lukas

% %% Change dissipation
% alpha=ssModel.A;
% %alpha=alpha*[0.8:0.05:1];  % comment for single simulation
% %alpha=[alpha,1];           % comment for single simulation
% 
% beta=ssModel.Bu;
% beta=beta*[0.8:0.05:1];     % comment for single simulation
% beta=[beta,1];              % comment for single simulation
% 
% for i = 1:1%length(alpha) %loop over alpha
% for itBeta = 1:length(beta)
% 
% % Reset constraints and objective
% con = [];
% obj = 0;
% 
% % New decision varaibles
% %s1 = sdpvar(6,N,'full');
% c = sdpvar(1, N,'full'); % CHF/kWh
% e = sdpvar(1, N,'full'); 
% v = sdpvar(1, N,'full'); 
% xb=sdpvar(1, N,'full');
% sb = sdpvar(1, N,'full'); 
% % Exercise specific parameters
% penal=1;
% 
% % Define constraints and objective for MPC-controller
% 
% for j = 1:N-1  
%     % System
%     obj = obj + c(j)*e(j)+  penal*S(:,j)'*S(:,j);   % Cost function
%     con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
%     con = [con, y(:,j) == C*x(:,j)];
%     con = [con, Hu*u(:,j) <= hu];                   % Input constraints
%     con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
%    
%     % Battery
%     con = [con, v(j)==e(j)-sum(u(:,j))];
%     con = [con, e(j)>=0];
%     con = [con, xb(j+1) == alpha(i)*xb(j)+beta(itBeta)*v(j)];
% 
% 
%     con = [con, -20 <= v(j) <= 20];  
%     
%     if j ~= 1
%         con = [con,  0 <=xb(j) <= 20];
%     end
% end
% 
% % Final constraints
%     % System
%     obj = obj + c(N)*e(N)+  penal*S(:,N)'*S(:,N);   % Cost function
%     con = [con, y(:,N) == C*x(:,N)];
%     con = [con, Hu*u(:,N) <= hu];                   % Input constraints
%     con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
%     % Battery
%     con = [con, v(N)==e(N)-sum(u(:,N))];
%     con = [con, e(N)>=0];
%     con = [con, 0 <= xb(N) <= 20];
%     con = [con, -20 <= v(N) <= 20];  
% 
%     % v power to battery 
%     % e power from grid
%     % u power to buildings
%     
% 
% ops = sdpsettings('verbose',1, 'solver', '+gurobi');
% controller = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);
% 
% [xt_bat, yt_bat, ut_bat, t_bat, et_bat, xbt_bat] = simBuildStorage(controller, T, @shiftPred, N,'fig/batt9');
% 
% axt_bat(:,:,i)=xt_bat;
% ayt_bat(:,:,i)=yt_bat;
% aut_bat(:,:,i)=ut_bat;
% at_bat(:,:,i)=t_bat;
% aet_bat(:,:,i)=et_bat;
% axbt_bat(:,:,i)= xbt_bat;
% 
% bxt_bat(:,:,itBeta)=xt_bat;
% byt_bat(:,:,itBeta)=yt_bat;
% but_bat(:,:,itBeta)=ut_bat;
% bt_bat(:,:,itBeta)=t_bat;
% bet_bat(:,:,itBeta)=et_bat;
% bxbt_bat(:,:,itBeta)= xbt_bat;
% 
% 
% %Total_bat(i)=refCost(1:T)/3*et_bat(1,:)';
% end
% end
% 
% %%
% if(length(alpha)>1)
%     % Plotting
%     figure
%     title('State of charge with respect to dissipation')
%     grid on
%     hold on
%     for i = 1:length(alpha)
%     plot(at_bat(:,:,i),axbt_bat(1,:,i),'LineWidth',1); 
%     hold on
%     end
%     legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');
% 
%     % Plotting
%     figure
%     title('Electricity from grid with respect to dissipation')
%     grid on
%     hold on
%     for i = 1:length(alpha)
%     plot(at_bat(:,:,i),aet_bat(1,:,i),'LineWidth',1);
%     hold on
%     aTotal_Cost_bat(i)=refCost(1:T)/3*aet_bat(1,:,i)';
%     end
%     legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');
% end
% 
% if(length(beta)>1)
%     % Plotting
%     figure
%     title('State of charge with respect to dissipation')
%     grid on
%     hold on
%     for i = 1:length(beta)
%     plot(bt_bat(1,:,i),bxbt_bat(1,:,i),'LineWidth',1); 
%     hold on
%     end
%     legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');
% 
%     % Plotting
%     figure
%     title('Electricity from grid with respect to dissipation')
%     grid on
%     hold on
%     for i = 1:length(beta)
%     plot(bt_bat(:,:,i),bet_bat(1,:,i),'LineWidth',1);
%     hold on
%     aTotal_Cost_bat(i)=refCost(1:T)/3*bet_bat(1,:,i)';
%     end
%     legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');
% end
% 
% 
% 
% %%
% % Plotting
% 
% figure
% title('State of charge with respect to varying storage capacity')
% grid on
% hold on
% for i = 1:length(Stor_cap)
% plot(bt_bat(:,:,i),bxbt_bat(1,:,i),'LineWidth',1); 
% hold on
% end
% legend('10 kWh','15 kWh','20 kWh','25 kWh','30 kWh','35 kWh')
% figure
% title('Electricity from grid with respect to varying storage capacity')
% grid on
% hold on
% for i = 1:length(Stor_cap)
% plot(bt_bat(:,:,i),bet_bat(1,:,i),'LineWidth',1);
% hold on
% bTotal_Cost_bat(i)=refCost(1:T)/3*bet_bat(1,:,i)';
% end
% legend('10 kWh','15 kWh','20 kWh','25 kWh','30 kWh','35 kWh')

% %% Carge dissipation
% alpha=ssModel.A;
% beta=ssModel.Bu;
% beta=beta*[0.5:0.05:1];
% 
% 
% for i = 1:length(beta)
% % Reset constraints and objective
% con = [];
% obj = 0;
% 
% % New decision varaibles
% %s1 = sdpvar(6,N,'full');
% c = sdpvar(1, N,'full'); % CHF/kWh
% e = sdpvar(1, N,'full'); 
% v = sdpvar(1, N,'full'); 
% xb=sdpvar(1, N,'full');
% sb = sdpvar(1, N,'full'); 
% % Exercise specific parameters
% penal=1;
% 
% 
% % Define constraints and objective for MPC-controller
% 
% for j = 1:N-1  
%     % System
%     obj = obj + c(j)*e(j)+  penal*S(:,j)'*S(:,j);   % Cost function
%     con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
%     con = [con, y(:,j) == C*x(:,j)];
%     con = [con, Hu*u(:,j) <= hu];                   % Input constraints
%     con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
%    
%     % Battery
%     con = [con, v(j)==e(j)-sum(u(:,j))];
%     con = [con, e(j)>=0];
%     con = [con, xb(j+1) == alpha*xb(j)+beta(i)*v(j)];
% 
% 
%     con = [con, -20 <= v(j) <= 20];  
%     
%     if j ~= 1
%         con = [con,  0 <=xb(j) <= 20];
%     end
% end
% 
% % Final constraints
%     % System
%     obj = obj + c(N)*e(N)+  penal*S(:,N)'*S(:,N);   % Cost function
%     con = [con, y(:,N) == C*x(:,N)];
%     con = [con, Hu*u(:,N) <= hu];                   % Input constraints
%     con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
%     % Battery
%     con = [con, v(N)==e(N)-sum(u(:,N))];
%     con = [con, e(N)>=0];
%     con = [con, 0 <= xb(N) <= 20];
%     con = [con, -20 <= v(N) <= 20];  
% 
%     % v power to battery 
%     % e power from grid
%     % u power to buildings
%    
%     
% ops = sdpsettings('verbose',1, 'solver', '+gurobi');
% controller = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);
% 
% [xt_bat, yt_bat, ut_bat, t_bat, et_bat, xbt_bat] = simBuildStorage(controller, T, @shiftPred, N,'hello');
% 
% axt_bat(:,:,i)=xt_bat;
% ayt_bat(:,:,i)=yt_bat;
% aut_bat(:,:,i)=ut_bat;
% at_bat(:,:,i)=t_bat;
% aet_bat(:,:,i)=et_bat;
% axbt_bat(:,:,i)= xbt_bat;
% 
% %Total_bat(i)=refCost(1:T)/3*et_bat(1,:)';
% end
% 
% % Plotting
% figure
% title('State of charge with respect to dissipation')
% grid on
% hold on
% for i = 1:length(beta)
% plot(at_bat(:,:,i),axbt_bat(1,:,i),'LineWidth',1); 
% hold on
% end
% legend('+50% dissipation','+45% dissipation','+40% dissipation','+35% dissipation','+30% dissipation','+25% dissipation','+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','+original dissipation');
% 
% % Plotting
% figure
% title('Electricity from grid with respect to dissipation2')
% grid on
% for i = 1:length(beta)
% plot(at_bat(:,:,i),aet_bat(1,:,i),'LineWidth',1);
% hold on
% end
% legend('+50% dissipation','+45% dissipation','+40% dissipation','+35% dissipation','+30% dissipation','+25% dissipation','+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','+original dissipation');

%% Exercise 5 Jannik
%% Change dissipation alpha
alpha=ssModel.A;
alpha=alpha*[0.8:0.05:1];
alpha=[alpha,1];

beta=ssModel.Bu;

for i = 1:length(alpha)
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
%s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
e = sdpvar(1, N,'full'); 
v = sdpvar(1, N,'full'); 
xb=sdpvar(1, N,'full');
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=1;


% Define constraints and objective for MPC-controller

for j = 1:N-1  
    % System
    obj = obj + c(j)*e(j)+  penal*S(:,j)'*S(:,j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
   
    % Battery
    con = [con, v(j)==e(j)-sum(u(:,j))];
    con = [con, e(j)>=0];
    con = [con, xb(j+1) == alpha(i)*xb(j)+beta*v(j)];


    con = [con, -20 <= v(j) <= 20];  
    
    if j ~= 1
        con = [con,  0 <=xb(j) <= 20];
    end
end

% Final constraints
    % System
    obj = obj + c(N)*e(N)+  penal*S(:,N)'*S(:,N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
    % Battery
    con = [con, v(N)==e(N)-sum(u(:,N))];
    con = [con, e(N)>=0];
    con = [con, 0 <= xb(N) <= 20];
    con = [con, -20 <= v(N) <= 20];  

    % v power to battery 
    % e power from grid
    % u power to buildings
    

ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);

[xt_bat, yt_bat, ut_bat, t_bat, et_bat, xbt_bat] = simBuildStorageJZ(controller, T, @shiftPred, N,'hello');

axt_bat(:,:,i)=xt_bat;
ayt_bat(:,:,i)=yt_bat;
aut_bat(:,:,i)=ut_bat;
at_bat(:,:,i)=t_bat;
aet_bat(:,:,i)=et_bat;
axbt_bat(:,:,i)= xbt_bat;

%Total_bat(i)=refCost(1:T)/3*et_bat(1,:)';
end

% Plotting
figure
title('State of charge with respect to dissipation')
grid on
hold on
for i = 1:length(alpha)
plot(at_bat(:,:,i),axbt_bat(1,:,i),'LineWidth',1); 
hold on
end
legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');

% Plotting
figure
title('Electricity from grid with respect to dissipation')
grid on
hold on
for i = 1:length(alpha)
plot(at_bat(:,:,i),aet_bat(1,:,i),'LineWidth',1);
hold on
aTotal_Cost_bat(i)=refCost(1:T)/3*aet_bat(1,:,i)';
end
legend('+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','original dissipation','no dissipation');

%% Change storage capacity

alpha=ssModel.A;
beta=ssModel.Bu;

Stor_cap=[10,15,20,25,30,35];


for i = 1:length(Stor_cap)
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
%s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
e = sdpvar(1, N,'full'); 
v = sdpvar(1, N,'full'); 
xb=sdpvar(1, N,'full');
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=1;


% Define constraints and objective for MPC-controller

for j = 1:N-1  
    % System
    obj = obj + c(j)*e(j)+  penal*S(:,j)'*S(:,j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
   
    % Battery
    con = [con, v(j)==e(j)-sum(u(:,j))];
    con = [con, e(j)>=0];
    con = [con, xb(j+1) == alpha*xb(j)+beta*v(j)];


    con = [con, -20 <= v(j) <= 20];  
    
    if j ~= 1
        con = [con,  0 <=xb(j) <= Stor_cap(i)];
    end
end

% Final constraints
    % System
    obj = obj + c(N)*e(N)+  penal*S(:,N)'*S(:,N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
    % Battery
    con = [con, v(N)==e(N)-sum(u(:,N))];
    con = [con, e(N)>=0];
    con = [con, 0 <= xb(N) <= 20];
    con = [con, -20 <= v(N) <= 20];  

    % v power to battery 
    % e power from grid
    % u power to buildings
    

ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);

[xt_bat, yt_bat, ut_bat, t_bat, et_bat, xbt_bat] = simBuildStorageJZ(controller, T, @shiftPred, N,'hello');

bxt_bat(:,:,i)=xt_bat;
byt_bat(:,:,i)=yt_bat;
but_bat(:,:,i)=ut_bat;
bt_bat(:,:,i)=t_bat;
bet_bat(:,:,i)=et_bat;
bxbt_bat(:,:,i)= xbt_bat;

end

% Plotting
figure
title('State of charge with respect to varying storage capacity')
grid on
hold on
for i = 1:length(Stor_cap)
plot(bt_bat(:,:,i),bxbt_bat(1,:,i),'LineWidth',1); 
hold on
end
legend('10 kWh','15 kWh','20 kWh','25 kWh','30 kWh','35 kWh')
figure
title('Electricity from grid with respect to varying storage capacity')
grid on
hold on
for i = 1:length(Stor_cap)
plot(bt_bat(:,:,i),bet_bat(1,:,i),'LineWidth',1);
hold on
bTotal_Cost_bat(i)=refCost(1:T)/3*bet_bat(1,:,i)';
end
legend('10 kWh','15 kWh','20 kWh','25 kWh','30 kWh','35 kWh')

%% Carge dissipation beta
alpha=ssModel.A;
beta=ssModel.Bu;
beta=beta*[0.5:0.05:1];


for i = 1:length(beta)
% Reset constraints and objective
con = [];
obj = 0;

% New decision varaibles
%s1 = sdpvar(6,N,'full');
c = sdpvar(1, N,'full'); % CHF/kWh
e = sdpvar(1, N,'full'); 
v = sdpvar(1, N,'full'); 
xb=sdpvar(1, N,'full');
sb = sdpvar(1, N,'full'); 
% Exercise specific parameters
penal=1;


% Define constraints and objective for MPC-controller

for j = 1:N-1  
    % System
    obj = obj + c(j)*e(j)+  penal*S(:,j)'*S(:,j);   % Cost function
    con = [con, x(:,j+1) == A*x(:,j) + Bu*u(:,j)+Bd*d(:,j)]; % System dynamics
    con = [con, y(:,j) == C*x(:,j)];
    con = [con, Hu*u(:,j) <= hu];                   % Input constraints
    con = [con, Hy*y(:,j) <= hy + S(:,j)+[sb(j);sb(j);sb(j);sb(j);sb(j);sb(j)]];% Output constraints
   
    % Battery
    con = [con, v(j)==e(j)-sum(u(:,j))];
    con = [con, e(j)>=0];
    con = [con, xb(j+1) == alpha*xb(j)+beta(i)*v(j)];


    con = [con, -20 <= v(j) <= 20];  
    
    if j ~= 1
        con = [con,  0 <=xb(j) <= 20];
    end
end

% Final constraints
    % System
    obj = obj + c(N)*e(N)+  penal*S(:,N)'*S(:,N);   % Cost function
    con = [con, y(:,N) == C*x(:,N)];
    con = [con, Hu*u(:,N) <= hu];                   % Input constraints
    con = [con, Hy*y(:,N) <= hy + S(:,N)+[sb(N);sb(N);sb(N);sb(N);sb(N);sb(N)]];% Output constraints
    % Battery
    con = [con, v(N)==e(N)-sum(u(:,N))];
    con = [con, e(N)>=0];
    con = [con, 0 <= xb(N) <= 20];
    con = [con, -20 <= v(N) <= 20];  

    % v power to battery 
    % e power from grid
    % u power to buildings
   
    
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c(:);sb(:)],[u;v;e]);

[xt_bat, yt_bat, ut_bat, t_bat, et_bat, xbt_bat] = simBuildStorageJZ(controller, T, @shiftPred, N,'hello');

axt_bat(:,:,i)=xt_bat;
ayt_bat(:,:,i)=yt_bat;
aut_bat(:,:,i)=ut_bat;
at_bat(:,:,i)=t_bat;
aet_bat(:,:,i)=et_bat;
axbt_bat(:,:,i)= xbt_bat;

%Total_bat(i)=refCost(1:T)/3*et_bat(1,:)';
end

% Plotting
figure
title('State of charge with respect to dissipation')
grid on
hold on
for i = 1:length(beta)
plot(at_bat(:,:,i),axbt_bat(1,:,i),'LineWidth',1); 
hold on
end
legend('+50% dissipation','+45% dissipation','+40% dissipation','+35% dissipation','+30% dissipation','+25% dissipation','+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','+original dissipation');

% Plotting
figure
title('Electricity from grid with respect to dissipation2')
grid on
for i = 1:length(beta)
plot(at_bat(:,:,i),aet_bat(1,:,i),'LineWidth',1);
hold on
end
legend('+50% dissipation','+45% dissipation','+40% dissipation','+35% dissipation','+30% dissipation','+25% dissipation','+20% dissipation','+15% dissipation','+10% dissipation','+5% dissipation','+original dissipation');






