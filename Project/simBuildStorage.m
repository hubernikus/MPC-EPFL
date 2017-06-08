%INPUTS:
%  controller - controller (optimizer instance) - the parameter order suggested for the controller is:
%           [x ; xb; d_pred(:) ; cp(:) ; sb(:)]            for economic MPC, where d_pred is the prediction of disturbance input d, cp is the electricity price vector, and sb is the night-setback offset, over the prediction horizon          
%  T - simulation time (in time-steps)
%  fhandle - function handle to the shiftPred function. The input to this
%               function is the time-step and the outputs of this function are the
%               predictions of disturbance d, electricity price cp, and the night-setback offset over the prediction horizon for the given time step. For
%               more details, check the documentation of this function.
%  N - Prediction Horizon of your MPC controller
%

%OUTPUTS:
% xt - state as a function of time
% yt - output as a function of time
% ut - input as a function of time
% t - time (time-steps)
% et - electricity input to the battery model as a function of time
% xbt - State of the battery storage as a function of time


function [ xt, yt, ut, t, et, xbt ] = simBuildStorage(controller, T, fhandle, N, figName)
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


x = x0red;
xb = 0;

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);


xt = zeros(nx,T);   % T is the length of simulation in samples
yt = zeros(ny,T);

ut = zeros(nu,T);
t = zeros(1,T);


et = zeros(1,T);
xbt = zeros(1,T);
vt = zeros(1,T);

cpt = zeros(1,T);
sbt = zeros(1,T);

% New Variables

% Variable Cost
Thours=T/3;
Tdays=Thours/24;
refCost =0.2*ones(1,length(refDist));
for i=0:floor(Tdays)  
refCost(i*24*3+(30:30+18))=0.04;
end

% Night Setback
setb=4*ones(1,length(refDist));
for i=0:floor(Tdays)  
setb(i*24*3+(24:24+30))=0;
end


for i = 1:T
[d_pred, cp, sb] = fhandle(i, N,refCost,setb);
[U, id] = controller{[x; xb; d_pred(:); cp(:); sb(:)]}; % this is the suggested form for the controller : you can change it provided buildSim.m is also accordingly changed

xt(:,i) = x;
ut(:,i) = U(1:nu,1);
et(:,i) = U(end,1);
vt(:,i) = U(end-1,1);
xbt(:,i) = xb;

cpt(:,i) = cp(1,1);
sbt(:,i) = sb(1,1);

yt(:,i) = C*x;
t(1,i) = i;

disp(['Iteration ' int2str(i)]);
yalmiperror(id);


x = A*x + Bu*ut(:,i) + Bd*d_pred(:,1);

xb = a*xb + b*[U(nu+1,1)];

end

%% Generating the Plots

% Converting time scale from time-step to hours
t = t./3;

figure('Position',[0 0 1200 500]);
% subplot(2,3,1)
plot(t, yt(1,:))
hold on; grid on
% plot(t, 26+sbt(1,:),'r')
% plot(t, 22-sbt(1,:),'r')
% legend('Zone-1 Temperature', 'Temperature Constraints')
xlabel('Hours');
ylabel('Temperature - Zone1 (C)');


% figure
% subplot(2,3,2)
plot(t, yt(2,:),'k')
% hold on
% plot(t, 26+sbt(1,:),'r')
% plot(t, 22-sbt(1,:),'r')
% legend('Zone-2 Temperature', 'Temperature Constraints')
xlabel('Hours');
ylabel('Temperature - Zone2 (C)');


% figure
% subplot(2,3,3)
plot(t, yt(3,:),'c')
hold on grid on;
plot(t, 26+sbt(1,:),'r')
plot(t, 22-sbt(1,:),'r')
% legend('Zone-3 Temperature', 'Temperature Constraints')
legend('Zone-1','Zone-2','Zone-3', 'Temperature Constraints')
xlabel('Hours');
ylabel('Temperature - Zone3 (C)');
xlim([t(1),t(end)])
print(sprintf('%s_temp',figName),'-dpng');

figure('Position',[0 0 1200 500]);
% subplot(2,3,4)
plot(t,ut(1,:))
hold on; grid on
% plot(t,10*cpt(1,:),'r')
% legend('Zone-1 Input', 'High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone1 (kW)');


% figure
% subplot(2,3,5)
plot(t,ut(2,:),'b'); grid on;
% hold on
% plot(t,10*cpt(1,:),'r')
% legend('Zone-2 Input', 'High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone2 (kW)');


% figure
% subplot(2,3,6)
plot(t,ut(3,:),'c')
hold on; grid on
plot(t,10*cpt(1,:),'r')
% legend('Zone-3 Input', 'High/Low Price Time')
legend('Zone-1','Zone-2','Zone-3','High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone3 (kW)');
xlim([t(1),t(end)])
print(sprintf('%s_powerZones',figName),'-dpng');

figure('visible','off');
subplot(2,1,1)
plot(t,xbt(1,:)); grid on
xlabel('Hours');
ylabel('Storage State');
xlim([t(1),t(end)])

% figure
subplot(2,1,2)
plot(t,et(1,:))
hold on; grid on;
plot(t,10*cpt(1,:),'r')
legend('Electrical Power Purchased','High/Low Price Time')
xlabel('Hours');
ylabel('Power purchased (kW)');
xlim([t(1),t(end)])
print(sprintf('%s_powerPurchase',figName),'-dpng');

figure('Position',[0 0 1200 800]);
plot(t,vt(1,:))
hold on; grid on;
plot(t,10*cpt(1,:),'r')
legend('Power to storage Purchased','High/Low Price Time')
xlabel('Hours');
ylabel('Power input to the Storage (kW)');
xlim([t(1),t(end)])
print(sprintf('%s_powerToStorage',figName), '-dpng');

end

