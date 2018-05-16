% loads log posterior and plots it
clear
file2load = 'data/logpost_tr100sd1part500sp50.mat';
load(file2load)
unnorm_post=exp(log_posterior);
samples=linspace(0,5,50);
true_h = 1;

% plotting parameters
lw=3;
fs=20;
plot(samples, unnorm_post, 'LineWidth', lw)
hold on
plot([true_h, true_h],[0, max(unnorm_post)], 'LineWidth',lw)
hold off
legend('posterior','true h')
title('posterior for stochastic nonlin')
ylabel('unnormalized posterior')
xlabel('h fit')
ax=gca; ax.FontSize=fs;
