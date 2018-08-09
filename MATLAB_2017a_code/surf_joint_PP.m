%% NL-L
clear
load('../data/joint_PP_ntrials_99999_noise_1.mat')
[X,Y]=meshgrid(gammas,hs);
%surf(X,Y,PP)
xlabel('gamma')
ylabel('h')
zlabel('predictive power')
% best gamma = 5.65
% best h = 0.45

[A,B]=max(PP);
hmax=zeros(size(B));
for i=1:length(B)
    hmax(i)=hs(B(i));
end
[C,D]=max(PP,[],2);
gmax=zeros(size(D));
for i=1:length(D)
    gmax(i)=gammas(D(i));
end
figure()
subplot(2,1,1)
plot(gammas,hmax)
xlabel('gamma')
ylabel('h_{max}')
subplot(2,1,2)
plot(hs,gmax)
xlabel('h')
ylabel('\gamma_{max}')

%% L-L
clear
load('../data/joint_PP_LL_ntrials_100_noise_1.mat')
[X,Y]=meshgrid(thetas_1,thetas_2);
%num_gammas=length(thetas_1);
surf(X,Y,PP)
xlabel('gamma')
ylabel('h')
zlabel('predictive power')

%% NL-NL
clear
load('../data/joint_PP_NLNL_ntrials_100_noise_1.mat')
[X,Y]=meshgrid(thetas_1,thetas_2);
%num_gammas=length(thetas_1);
surf(X,Y,PP)
xlabel('gamma')
ylabel('h')
zlabel('predictive power')