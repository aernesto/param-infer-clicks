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

[C,D]=max(PP');

figure()
subplot(2,1,1)
plot(gammas,A)
xlabel('gamma')
ylabel('h - best choice match')
subplot(2,1,2)
plot(hs,C)
xlabel('h')
ylabel('gamma - best choice match')