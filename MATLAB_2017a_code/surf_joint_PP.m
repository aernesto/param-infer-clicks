clear
load('../data/joint_PP_ntrials_99999_noise_1.mat')
[X,Y]=meshgrid(gammas,hs);
surf(X,Y,PP)
xlabel('gamma')
ylabel('h')
zlabel('predictive power')
% best gamma = 5.65
% best h = 0.45