% produces plot from data generated by 'explore_PP_L.m'
clear
tot_num_trials=100000; 
nsd=2;
gcode='00.11.5';
load(['../data/explore_PP_L_tr',num2str(tot_num_trials),'noise',...
    num2str(nsd),'gcode',gcode,'.mat'])

figure();
plot(gammas,PP)
title(['L ; noise=',num2str(nsd),' ; NumTrials=',...
    num2str(tot_num_trials)])
ylabel('Pred. Pow.')
xlabel('gamma')
ax=gca; ax.FontSize=16;