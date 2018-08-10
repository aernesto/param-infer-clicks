%% NL-L
clear
fs=10; % font size
lw=2; % linewidth
p=load('../data/joint_PP_ntrials_99999_noise_1.mat');
m=load(['../data/mse_nonlin_fig4_iteration2_',num2str(500),'trials.mat']);
pp_NLL=mapmode2pp(m.modes_nonlinlin(m.modes_nonlinlin<=2.5), p.hs, p.PP(:,68));
figure()
boxplot(pp_NLL)

hold on
ns=[68];

for ii=1
%     ax=subplot(2,1,ii);

    n=ns(ii);
%     plot(p.hs,p.PP(:,n),'LineWidth',lw)
    ax=gca;
    [MX,mx]=max(p.PP(:,n));
%     hold on
    if ii==1
        c1='--b'; c2='-b';
    else
        c1='--r'; c2='-r';
    end
    plot([.5,1.5],[MX,MX],'LineWidth',lw-1)
%     hold off
    title(['100K trials - \gamma = ',num2str(p.gammas(n))])
    xlabel('h')
    ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
%     ax.FontSize=fs;
end

hold off


%% L-NL
clear
fs=10; % font size
lw=2; % linewidth
load('../data/joint_PP_ntrials_99999_noise_1.mat')

figure()
hold on


%     ax=subplot(2,1,ii);

n=11;
plot(gammas,PP(n,:),'LineWidth',lw)
ax=gca;
[~,mx]=max(PP(n,:));
%     hold on

c1='--b'; c2='-b';

plot([gammas(mx),gammas(mx)],[ax.YLim(1),ax.YLim(2)],c2,...
    'LineWidth',lw-1)
%     hold off
title(['100K trials - h = ',num2str(hs(n))])
xlabel('\gamma')
ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
%     ax.FontSize=fs;

hold off

% L-L
clear
fs=10; % font size
lw=2; % linewidth

load('../data/joint_PP_LL_ntrials_999990_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_10000_noise_1.mat')

[X,Y]=meshgrid(thetas_1,thetas_2);
num_g=length(thetas_1);

for i=1:num_g
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end

figure()
hold on
ns=[68,69];
for ii=1:2
%     ax=subplot(2,1,ii);

    n=ns(ii);
    plot(thetas_1,PP(n,:),'LineWidth',lw)
    ax=gca;
    [~,mx]=max(PP(n,:));
%     hold on
    if ii==1
        c1='--b'; c2='-b';
    else
        c1='--r'; c2='-r';
    end
    plot([thetas_1(n),thetas_1(n)],[ax.YLim(1),ax.YLim(2)],c1,...
        'LineWidth',lw)
    plot([thetas_1(mx),thetas_1(mx)],[ax.YLim(1),ax.YLim(2)],c2,...
        'LineWidth',lw-1)
%     hold off
    title(['1M trials - ref \theta_1 = ',num2str(thetas_1(n))])
    xlabel('\theta_2')
    ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
%     ax.FontSize=fs;
end
hold off

% NL-NL
clear
fs=10; % font size
lw=2; % linewidth

load('../data/joint_PP_NLNL_ntrials_999990_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_10000_noise_1.mat')

[X,Y]=meshgrid(thetas_1,thetas_2);
num_g=length(thetas_1);

for i=1:num_g
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end

figure()
hold on
ns=[11];
for ii=1
%     ax=subplot(2,1,ii);

    n=ns(ii);
    plot(thetas_1,PP(n,:),'LineWidth',lw)
    ax=gca;
    [~,mx]=max(PP(n,:));
%     hold on
    if ii==1
        c1='--b'; c2='-b';
    else
        c1='--r'; c2='-r';
    end
    plot([thetas_1(n),thetas_1(n)],[ax.YLim(1),ax.YLim(2)],c1,...
        'LineWidth',lw)
    plot([thetas_1(mx),thetas_1(mx)],[ax.YLim(1),ax.YLim(2)],c2,...
        'LineWidth',lw-1)
%     hold off
    title(['NLNL- 1M trials - ref \theta_1 = ',num2str(thetas_1(n))])
    xlabel('\theta_2')
    ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
%     ax.FontSize=fs;
end
hold off

% NEW WHISKERS
% plot whiskers of accuracy recovered from modes of posteriors for
% fits of stochastic models
clear
for TN=[100,500]  % loop over block size (in trial numbers)
    % nonlin
    % load modes
    load(['../data/mse_nonlin_fig4_iteration2_',num2str(TN),'trials.mat'])
    load(['../data/mse_lin_fig4_iteration2_',num2str(TN),'trials.mat'])
    % load accuracy values (obtained from my .h5 DBs)
    %load('../data/linacc.csv','-ascii')  % loads linacc
    load('../data/linacc_S3lr2.csv','-ascii')  % loads linacc for db S3lr2
    linacc=linacc_S3lr2;
    
    load('../data/nonlinacc.csv','-ascii')  % loads nonlinacc
    accuracy_nonlinlin = mapmode2pp(modes_nonlinlin, linspace(0,10,1000),...
        nonlinacc');
    accuracy_nonlinnonlin = mapmode2pp(modes_nonlinnonlin,...
        linspace(0,10,1000), nonlinacc');
    accuracy_linlin = mapmode2pp(modes_linlin, linspace(0,40,10000),...
        linacc');
    accuracy_linnonlin = mapmode2pp(modes_linnonlin,...
        linspace(0,40,10000), linacc');
    %g=['NLL','NLNL','LL','LNL'];%[ones(500,1),2*ones(500,1),3*ones(500,1),4*ones(500,1)];
    lw=4;
    ms=8;
    fs=20;
    boxplot([accuracy_nonlinlin;...
        accuracy_nonlinnonlin;...
        accuracy_linlin;...
        accuracy_linnonlin]')
    set(findobj(gca,'type','line'),'linew',lw)
    set(gca,'linew',lw/2)
%     ylim([.85,.9])
    ylabel('PP')
    ax=gca;ax.FontSize=fs;
    saveas(gcf, ['whisker_PP_',num2str(TN),'.png'])
end

