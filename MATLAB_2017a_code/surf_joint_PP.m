%% NL-L
clear
fs=20; % font size
lw=3; % linewidth

load('../data/joint_PP_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LNL_ntrials_999990_noise_2.mat')
%load('../data/joint_PP_LNL_ntrials_999990_noise_0.1.mat')
%gammas=thetas_1;hs=thetas_2;
PP=PP';
[X,Y]=meshgrid(gammas,hs);
[A,B]=max(PP');
hmax=zeros(size(B));
for i=1:length(B)
    hmax(i)=hs(B(i));
end
[C,D]=max(PP',[],2);
gmax=zeros(size(D));
for i=1:length(D)
    gmax(i)=gammas(D(i));
end
figure()
ax1=subplot(2,1,1);
pcolor(ax1,X,Y,PP')
colormap('copper')
hold on
plot(ax1,gammas,hmax,'LineWidth',lw,'Color','red')
plot(ax1,hs,gmax,'LineWidth',lw,'Color','red')
hold off
title('Predictive Power NL-L')
xlabel('\gamma')
ylabel('h')
hold off
ax1.FontSize=fs;

%figure()
ax2=subplot(2,1,2);
surf(X,Y,PP')
title('noise = 0.1')
xlabel('\gamma')
ylabel('h')
zlabel('predictive power')
ax2.FontSize=fs;
% best gamma = 5.65
% best h = 0.45


% figure()
% subplot(2,1,1)
% plot(gammas,hmax)
% xlabel('gamma')
% ylabel('h_{max}')
% subplot(2,1,2)
% plot(hs,gmax)
% xlabel('h')
% ylabel('\gamma_{max}')

%% L-L
clear
fs=20; % font size
lw=3; % linewidth

%load('../data/joint_PP_LL_ntrials_999990_noise_0.1.mat')
load('../data/joint_PP_LL_ntrials_999990_noise_2.mat')
%load('../data/joint_PP_LL_ntrials_999990_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_LL_ntrials_10000_noise_1.mat')

[X,Y]=meshgrid(thetas_1,thetas_2);
num_g=length(thetas_1);

for i=1:num_g
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end


[A,B]=max(PP);
th1max=zeros(size(B));
for i=1:length(B)
    th1max(i)=thetas_1(B(i));
end

% figure()
% subplot(2,1,1)
% hold on
% plot(ax1,thetas_1,th1max,'LineWidth',lw)
% plot(ax1,th1max,thetas_1,'LineWidth',lw)
% plot(thetas_1,thetas_1,'--r','LineWidth',lw-1)
% hold off
% title('max PP - LL - 1M trials')
% xlabel('\theta_2')
% ylabel('\theta_1^{max}')
% ax=gca; ax.FontSize=fs;
% 
% subplot(2,1,2)
% plot(thetas_1,A)
% ylabel('max PP')



figure()
ax1=subplot(2,1,1);
pcolor(ax1,X,Y,PP)
colormap('copper')
hold on
plot(ax1,thetas_1,th1max,'LineWidth',lw,'Color','red')
plot(ax1,th1max,thetas_1,'LineWidth',lw,'Color','red')
plot(thetas_1,thetas_1,'--r','LineWidth',lw-1)
hold off
title('Predictive Power L-L')
xlabel('\gamma')
ylabel('\gamma')
hold off
ax1.FontSize=fs;

%figure()
ax2=subplot(2,1,2);
surf(X,Y,PP)
title('noise = 0.1')
xlabel('\gamma')
ylabel('\gamma')
zlabel('predictive power')
ax2.FontSize=fs;



% figure()
% for ii=1:3
%     ax=subplot(3,1,ii);
%     n=ii*floor(num_g/3);
%     plot(thetas_1,PP(n,:),'LineWidth',lw)
%     [~,mx]=max(PP(n,:));
%     hold on
%     plot([thetas_1(n),thetas_1(n)],[ax.YLim(1),ax.YLim(2)],'--r',...
%         'LineWidth',lw)
%     plot([thetas_1(mx),thetas_1(mx)],[ax.YLim(1),ax.YLim(2)],'-k',...
%         'LineWidth',lw-1)
%     hold off
%     title(['1M trials - ref \theta_1 = ',num2str(thetas_1(n))])
%     xlabel('\theta_2')
%     ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
% end


%% NL-NL
clear
fs=20; % font size
lw=3; % linewidth


%load('../data/joint_PP_NLNL_ntrials_10000_noise_1.mat')
%load('../data/joint_PP_NLNL_ntrials_99999_noise_1.mat')
%load('../data/joint_PP_NLNL_ntrials_999990_noise_1.mat')
load('../data/joint_PP_NLNL_ntrials_999990_noise_2.mat')
%load('../data/joint_PP_NLNL_ntrials_999990_noise_0.1.mat')


[X,Y]=meshgrid(thetas_1,thetas_2);
num_h=length(thetas_1);
for i=1:num_h
    for j=1:i-1
        PP(i,j)=PP(j,i);
    end
end

[A,B]=max(PP);
th1max=zeros(size(B));
for i=1:length(B)
    th1max(i)=thetas_1(B(i));
end

% figure()
% subplot(2,1,1)
% plot(thetas_1,h1max,'LineWidth',lw)
% hold on
% plot(thetas_1,thetas_1,'--r','LineWidth',lw-1)
% hold off
% 
% title('max PP - NLNL')
% xlabel('\theta_2')
% ylabel('\theta_1^{max}')
% ax=gca; ax.FontSize=fs;
% 
% subplot(2,1,2)
% plot(thetas_1,A)
% ylabel('max PP')


figure()
ax1=subplot(2,1,1);
pcolor(ax1,X,Y,PP)
colormap('copper')
hold on
plot(ax1,thetas_1,th1max,'LineWidth',lw,'Color','red')
plot(ax1,th1max,thetas_1,'LineWidth',lw,'Color','red')
plot(thetas_1,thetas_1,'--r','LineWidth',lw-1)
hold off
title('Predictive Power NL-NL')
xlabel('\gamma')
ylabel('\gamma')
hold off
ax1.FontSize=fs;

%figure()
ax2=subplot(2,1,2);
surf(X,Y,PP)
title('noise = 0.1')
xlabel('\gamma')
ylabel('\gamma')
zlabel('predictive power')
ax2.FontSize=fs;



% figure()
% for ii=1:3
%     ax=subplot(3,1,ii);
%     n=ii*floor(num_h/3);
%     plot(thetas_1,PP(n,:),'LineWidth',lw)
%     [~,mx]=max(PP(n,:));
%     hold on
%     plot([thetas_1(n),thetas_1(n)],[ax.YLim(1),ax.YLim(2)],'--r',...
%         'LineWidth',lw)
%     plot([thetas_1(mx),thetas_1(mx)],[ax.YLim(1),ax.YLim(2)],'-k',...
%         'LineWidth',lw-1)
%     hold off
%     title(['ref \theta_1 = ',num2str(thetas_1(n))])
%     xlabel('\theta_2')
%     ylabel('PP')
%     legend('PP','\theta_1=\theta_2','max')
%     ax.FontSize=fs;
% end