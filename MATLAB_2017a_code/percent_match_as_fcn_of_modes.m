% plot % match as fcn of mode, for each model pair

clear

for block_size=100:100:500
    bs_str=num2str(block_size);
    load(['../data/mse_nonlin_fig4_iteration2_',bs_str,'trials.mat'])
    load(['../data/mse_lin_fig4_iteration2_',bs_str,'trials.mat'])
    modes=[modes_linlin;modes_linnonlin;modes_nonlinnonlin;modes_nonlinlin]';
    
    
    percent_match=load(['../data/choice_match_',bs_str,'_4.mat']);
    
    lw=4;
    ms=8;
    fs=20;
    
    % compute correlation coefs
    LL_corr = corrcoef(percent_match.match(:,1),modes(:,1));
    LL_corr = LL_corr(2,1);
    
    LNL_corr = corrcoef(percent_match.match(:,2),modes(:,2));
    LNL_corr = LNL_corr(2,1);
    
    NLNL_corr = corrcoef(percent_match.match(:,3),modes(:,3));
    NLNL_corr = NLNL_corr(2,1);
    
    NLL_corr = corrcoef(percent_match.match(:,4),modes(:,4));
    NLL_corr = NLL_corr(2,1);
    
    
    plot(modes,percent_match.match,'o','MarkerSize',8)
    legend({['L-L R=',num2str(LL_corr)],...
        ['L-NL R=',num2str(LNL_corr)],...
        ['NL-NL R=',num2str(NLNL_corr)],...
        ['NL-L R=',num2str(NLL_corr)]},'Location','southoutside')
    ylim([0.78,.95])
    xlim([0,20])
    ax=gca;
    ax.FontSize=fs;
    saveas(gcf, ['../figs/match_fcn_modes_2',bs_str,'.pdf'])
end