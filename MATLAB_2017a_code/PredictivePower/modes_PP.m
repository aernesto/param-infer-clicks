%% NL-L
clear
fs=10; % font size
lw=2; % linewidth
p_NLL=load('../data/joint_PP_ntrials_99999_noise_1.mat');
m500=load(['../data/mse_nonlin_fig4_iteration2_',num2str(500),'trials.mat']);
m100=load(['../data/mse_nonlin_fig4_iteration2_',num2str(100),'trials.mat']);

fprintf('NLL 500 %d modes > 2.5\n',sum(m500.modes_nonlinlin>2.5))
fprintf('NLL 100 %d modes > 2.5\n',sum(m100.modes_nonlinlin>2.5))
pp_NLL_500=mapmode2pp(m500.modes_nonlinlin(m500.modes_nonlinlin<=2.5), p_NLL.hs, p_NLL.PP(:,68));
pp_NLL_100=mapmode2pp(m100.modes_nonlinlin(m100.modes_nonlinlin<=2.5), p_NLL.hs, p_NLL.PP(:,68));

figure()
n=68; nt=[100,500];
for ii=1:2
    ax=subplot(1,2,ii);
    if ii==1; boxplot(pp_NLL_100); else boxplot(pp_NLL_500); end
    ax=gca;
    [MX,mx]=max(p_NLL.PP(:,n));
    hold on
    c1='--b'; c2='-b';
    c1='--r'; c2='-r';
    plot([.5,1.5],[MX,MX],'LineWidth',lw-1)
    hold off
    title(['NL-L ',num2str(nt(ii)),' trials'])
    ylabel('PP')
    ylim([.84,.885])
%     legend('PP','\theta_1=\theta_2','max')
    ax.FontSize=fs;
end

%% L-NL
clear
fs=10; % font size
lw=2; % linewidth
p_LNL=load('../data/joint_PP_ntrials_99999_noise_1.mat');
m500=load(['../data/mse_lin_fig4_iteration2_',num2str(500),'trials.mat']);
m100=load(['../data/mse_lin_fig4_iteration2_',num2str(100),'trials.mat']);

fprintf('LNL 500 %d modes > 10\n',sum(m500.modes_linnonlin>10))
fprintf('LNL 100 %d modes > 10\n',sum(m100.modes_linnonlin>10))
pp_LNL_500=mapmode2pp(m500.modes_linnonlin(m500.modes_linnonlin<=10), p_LNL.gammas, p_LNL.PP(11,:));
pp_LNL_100=mapmode2pp(m100.modes_linnonlin(m100.modes_linnonlin<=10), p_LNL.gammas, p_LNL.PP(11,:));

figure()
n=11; nt=[100,500];
for ii=1:2
    ax=subplot(1,2,ii);
    if ii==1; boxplot(pp_LNL_100); else boxplot(pp_LNL_500); end
    ax=gca;
    [MX,mx]=max(p_LNL.PP(n,:));
    hold on
    c1='--b'; c2='-b';
    c1='--r'; c2='-r';
    plot([.5,1.5],[MX,MX],'LineWidth',lw-1)
    hold off
    title(['L-NL ',num2str(nt(ii)),' trials'])
    ylabel('PP')
    ylim([.84,.885])
%     legend('PP','\theta_1=\theta_2','max')
    ax.FontSize=fs;
end

%% L-L
clear
fs=10; % font size
lw=2; % linewidth

load('../data/joint_PP_LL_ntrials_999990_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_10000_noise_1.mat')
m500=load(['../data/mse_lin_fig4_iteration2_',num2str(500),'trials.mat']);
m100=load(['../data/mse_lin_fig4_iteration2_',num2str(100),'trials.mat']);

fprintf('LL 500 %d modes > 10\n',sum(m500.modes_linlin>10))
fprintf('LL 100 %d modes > 10\n',sum(m100.modes_linlin>10))

num_h=length(thetas_1);

for i=1:num_h
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end
hs=thetas_1;
pp_NLNL_500=mapmode2pp(m500.modes_linlin(m500.modes_linlin<=10), hs, PP(68,:));
pp_NLNL_100=mapmode2pp(m100.modes_linlin(m100.modes_linlin<=10), hs, PP(68,:));

figure()
n=68; nt=[100,500];
for ii=1:2
    ax=subplot(1,2,ii);
    if ii==1; boxplot(pp_NLNL_100); else boxplot(pp_NLNL_500); end
    ax=gca;
    [MX,mx]=max(PP(n,:));
    hold on
    c1='--b'; c2='-b';
    c1='--r'; c2='-r';
    plot([.5,1.5],[MX,MX],'LineWidth',lw-1)
    hold off
    title(['L-L ',num2str(nt(ii)),' trials'])
    ylabel('PP')
    ylim([.84,.885])
%     legend('PP','\theta_1=\theta_2','max')
    ax.FontSize=fs;
end



%% NL-NL
clear
fs=10; % font size
lw=2; % linewidth

load('../data/joint_PP_NLNL_ntrials_999990_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_10000_noise_1.mat')
m500=load(['../data/mse_nonlin_fig4_iteration2_',num2str(500),'trials.mat']);
m100=load(['../data/mse_nonlin_fig4_iteration2_',num2str(100),'trials.mat']);

fprintf('NLNL 500 %d modes > 2.5\n',sum(m500.modes_nonlinnonlin>2.5))
fprintf('NLNL 100 %d modes > 2.5\n',sum(m100.modes_nonlinnonlin>2.5))

num_h=length(thetas_1);

for i=1:num_h
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end
hs=thetas_1;
pp_NLNL_500=mapmode2pp(m500.modes_nonlinnonlin(m500.modes_nonlinnonlin<=2.5), hs, PP(11,:));
pp_NLNL_100=mapmode2pp(m100.modes_nonlinnonlin(m100.modes_nonlinnonlin<=2.5), hs, PP(11,:));

figure()
n=11; nt=[100,500];
for ii=1:2
    ax=subplot(1,2,ii);
    if ii==1; boxplot(pp_NLNL_100); else boxplot(pp_NLNL_500); end
    ax=gca;
    [MX,mx]=max(PP(n,:));
    hold on
    c1='--b'; c2='-b';
    c1='--r'; c2='-r';
    plot([.5,1.5],[MX,MX],'LineWidth',lw-1)
    hold off
    title(['NL-NL ',num2str(nt(ii)),' trials'])
    ylabel('PP')
    ylim([.84,.885])
%     legend('PP','\theta_1=\theta_2','max')
    ax.FontSize=fs;
end

