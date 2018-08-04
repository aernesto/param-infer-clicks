% explore predictive power
clear
tic
nbins=1000;  % nb of bins to bin q (to estimate its prob. density)
db='../data/validation2.h5';  % db from which trials should be taken
tot_num_trials=100000; 
nsd=2;

step=.1;
gammas = 0:step:1.5; num_gammas=length(gammas);
gcode=[num2str(gammas(1)),num2str(step),num2str(gammas(end))];
PP=zeros(size(gammas));
for g_idx=1:num_gammas
    PP(g_idx)=PredictivePower(nbins,db,tot_num_trials,gammas(g_idx),nsd);
end

% figure();
% plot(gammas,PP)
% title(['L model ; noise=',num2str(nsd),' ; num_trials=',...
%     num2str(tot_num_trials)])
% ylabel('Pred. Pow.')
% xlabel('gamma')
% ax=gca; ax.FontSize=20;
toc
save(['../data/explore_PP_L_tr',num2str(tot_num_trials),'noise',...
    num2str(nsd),'gcode',gcode,'.mat'])