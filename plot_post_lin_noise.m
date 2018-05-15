clear
dbname='data/S3lr2h1T2tr5Ksp10K.h5';
gname='/lr2hr14h1T2';
gammas=(0:.05:20)';
noises=linspace(.1,3,4);
legend_s=cell(1,length(noises)+1);
base_string = 'noise ';
trialcounts=[50,200,500];
for j=1:length(trialcounts)
    ntrials=trialcounts(j);
    probs=zeros(length(noises),length(gammas));
    i=0;
    for std = noises
        i=i+1;
        p = lhd_lin_mul_tr_gauss_clicks(dbname,...
            gname, ntrials, std, gammas);
        probs(i,:)=exp(p')/max(exp(p));%p'/max(p);
        legend_s{i}=[base_string, num2str(std)];
    end
    figure(j)
    probs(isnan(probs))=0;
    plot(repmat(gammas,1,length(noises)),probs','LineWidth',3)
    hold on
    true_gamma = h5readatt(dbname, [gname,'/decision_lin'], 'best_gamma');
    legend_s{end}=['true gamma ',num2str(true_gamma)];
    plot([true_gamma, true_gamma],ylim, 'k','LineWidth',3)
    hold off
    title(['stoch. lin. 2 lin. fit - ',num2str(ntrials),' trials'])
    xlabel('gamma')
    ylabel('scaled posterior')
    legend(legend_s)
    ax=gca; ax.FontSize=20;
end