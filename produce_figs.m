function produce_figs(dbname, dsetname, numtrials, init_g, last_g, numsamples, init_h, last_h, savelocation)
    % computes empirical probabilities of each sample being the correct one
    dsetname_decision_nonlin = [dsetname,'/decision_nonlin'];
    dsetname_decision_lin = [dsetname,'/decision_lin'];
    dsetname_info = [dsetname,'/trial_info'];
    linear_decisions = h5read(dbname, dsetname_decision_lin, [2,1], [numsamples, numtrials]);
    nonlinear_decisions = h5read(dbname, dsetname_decision_nonlin, [2,1], [numsamples, numtrials]);
    reference_decision_linear = h5read(dbname, dsetname_decision_lin, [1,1], [1, numtrials]);
    reference_decision_nonlinear = h5read(dbname, dsetname_decision_nonlin, [1,1], [1, numtrials]);
    true_param_lin = h5readatt(dbname,dsetname_decision_lin,'best_gamma');
    true_param_nonlin = h5readatt(dbname, dsetname_info, 'h');
    lin_param_values = linspace(init_g,last_g,numsamples);
    nonlin_param_values = linspace(init_h, last_h, numsamples)
    % linlin
    bool_lin_lin = (linear_decisions == reference_decision_linear);
    proba_linear = mean(bool_lin_lin, 2);
    figure('Visible', 'off')
    plot(lin_param_values,proba_linear, 'LineWidth', 3)
    hold on
    plot([true_param_lin true_param_lin],[0 1],'red','LineWidth', 2)
    hold off
    title('Lin fit to lin')
    xlabel('Disc. param. gamma')
    ylabel('% compatible')
    legend('proportion trials','true gamma (best)','Location','southeast')
    ax=gca;
    ax.FontSize=20;
    saveas(gcf,[savelocation,'fig_',dsetname(2:end),'_tr',...
                num2str(numtrials),'sp',num2str(numsamples),'_lin2lin.png'])

    % nonlinnonlin
    bool_nonlin_nonlin = (nonlinear_decisions == reference_decision_nonlinear);
    proba_nonlinear = mean(bool_nonlin_nonlin, 2);
    figure('Visible', 'off')
    plot(nonlin_param_values,proba_nonlinear, 'LineWidth', 3)
    hold on
    plot([true_param_nonlin true_param_nonlin],[0 1],'red','LineWidth', 2)
    hold off
    title('Nonlin fit to nonlin')
    xlabel('Discount. param. h')
    ylabel('% compatible')
    legend('proportion trials','true hazard rate','Location','southeast')
    ax=gca;
    ax.FontSize=20;
    saveas(gcf, [savelocation, 'fig_',dsetname(2:end),'_tr',...
                num2str(numtrials),'sp',num2str(numsamples), '_nonlin2nonlin.png'])
    % linnonlin - fit linear model to nonlinear data
    bool_lin_nonlin = (linear_decisions == reference_decision_nonlinear);
    proba_linnonlinear = mean(bool_lin_nonlin, 2);

    figure('Visible', 'off')
    plot(lin_param_values,proba_linnonlinear, 'LineWidth', 3)
    hold on
    plot([true_param_lin true_param_lin],[0 1],'red','LineWidth', 2)
    hold off
    title('Lin fit to nonlin')
    xlabel('Disc. param. gamma')
    ylabel('% compatible')
    legend('proportion trials','best gamma','Location','southeast')
    ax=gca;
    ax.FontSize=20;
    saveas(gcf,[savelocation,'fig_',dsetname(2:end),'_tr',...
                num2str(numtrials),'sp',num2str(numsamples),'_lin2nonlin.png'])
    % nonlinlin - fit nonlinear model to linear data
    bool_nonlin_lin = (nonlinear_decisions == reference_decision_linear);
    proba_nonlinlinear = mean(bool_nonlin_lin, 2);
    figure('Visible', 'off')
    plot(nonlin_param_values,proba_nonlinlinear,'LineWidth', 3)
    hold on
    plot([true_param_nonlin true_param_nonlin],[0 1],'red','LineWidth', 2)
    hold off
    title('Nonlin fit to lin')
    xlabel('Discount. param. h')
    ylabel('% compatible')
    legend('proportion trials','true h','Location','southeast')
    ax=gca;
    ax.FontSize=20;
    saveas(gcf,[savelocation, 'fig_',dsetname(2:end),'_tr',...
                num2str(numtrials),'sp',num2str(numsamples),'_nonlin2lin.png'])
end
