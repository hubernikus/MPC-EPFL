%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                       
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all;

clear all;
%%

xPos = [1 2 3 4 5];
cost = [566, 456, 389,251,232];
energ = [2.83,2.36,2.60,2.06,2.60];

figure('Position',[0 0 700 400])
yyaxis left; grid on;
plot(xPos,cost,'-o')
yyaxis right;
plot(xPos, energ,'-x')
legend('Cost [USD]','Energy Demand [kWh]')
xticks(xPos);
xticklabels({'Simple MPC','Economic','Variable Cost','Nighttime Setback','Battery Storage'})

print('fig/evolutionEnergyCost','-dpng')