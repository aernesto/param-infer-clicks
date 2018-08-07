% plots whiskers to illustrate spread of modes of posteriors
clear

lw=3;
ms=8;
fs=20;

true_h=1;
true_g= 6.7457;
nsd=1; 
for block_size=[1,20:20:80]
    load(['../data/fit_lin_stoch_',num2str(nsd),'_trials_',num2str(block_size)])
    %load(['../data/mse_nonlin_fig4_iteration2_',num2str(block_size),'trials.mat'])
    %load(['../data/mse_lin_fig4_iteration2_',num2str(block_size),'trials.mat'])
    modes=[modes_linlin;modes_linnonlin]';%modes_nonlinnonlin;modes_nonlinlin]';
    
    boxplot(modes)
    set(findobj(gca,'type','line'),'linew',lw)
    set(gca,'linew',lw/2)
    
    hold on
    ax=gca;
    plot([ax.XLim(1), ax.XLim(2)],[true_g,true_g],'LineWidth',lw/1.5)
%    plot([ax.XLim(1),ax.XLim(2)],[true_h,true_h],'LineWidth',lw/1.5)
    hold off
    
    %ax.XTickLabels=['L-L','L-NL','NL-NL','NL-L'];
    ylim([0,19])
    ax.FontSize=fs;
    saveas(gcf, ['../figs/whisker_modes_L_',num2str(block_size),'.pdf'])
end