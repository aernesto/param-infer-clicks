% plot MSE for nonlin-to-nonlin model fitting - using data from srvr
lw=3;
noise=1.5;
trial_nb = [50:50:500,600];
M=length(trial_nb);
filenames = {'data/mses500runs50trials.mat',...
    'data/mses500runs100trials.mat',...
    'data/mses500runs150trials.mat',...
    'data/mses400runs200trials.mat',...
    'data/mses400runs250trials.mat',...
    'data/mses400runs300trials.mat',...
    'data/mses400runs350trials.mat',...
    'data/mses250runs400trials.mat',...
    'data/mses250runs450trials.mat',...
    'data/mses250runs500trials.mat',...
    'data/mses250runs600trials.mat'};
mse_vec = zeros(1,M);
for i = 1:M
    load(filenames{i})
    mse_vec(i)=mses;
end
plot(trial_nb, mse_vec,...
    'LineWidth',lw,... 
    'LineStyle','-o-')
xlabel('trial nb')
ylabel('MSE')
title(['nonlinnonlin noise=',num2str(noise)])