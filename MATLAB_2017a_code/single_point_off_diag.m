% reproducing the off-diag result for a single point on S3lr5h1T2tr...
clear
tic
rng('shuffle')
g1 = 6;
g2 = 5;

% load clicks
filename='../data/S3lr5h1T2tr10000sp1000.h5';
num_trials=10000;
[tr,envt]=get_trials(filename, num_trials);
kappa=log(20/5);
noise_sd=1; y0=0;
num_runs=1000;
runs=zeros(num_runs,2);
for run=1:num_runs
    match1=0;match2=0;
    for trn = 1:num_trials
        [lst, rst]=tr{:,trn};
        dec_g1_1=gauss_noise_lin_decide(lst, rst,g1, kappa, noise_sd, y0);
        dec_g1_2=gauss_noise_lin_decide(lst, rst,g1, kappa, noise_sd, y0);
        dec_g2_1=gauss_noise_lin_decide(lst, rst, g2, kappa, noise_sd, y0);
        dec_g2_2=gauss_noise_lin_decide(lst, rst, g1, kappa, noise_sd, y0);
        if dec_g1_1==0; dec_g1_1=sign(rand-0.5); end
        if dec_g1_2==0; dec_g1_2=sign(rand-0.5); end
        if dec_g2_1==0; dec_g2_1=sign(rand-0.5); end
        if dec_g2_2==0; dec_g2_2=sign(rand-0.5); end
        if dec_g1_1==dec_g1_2; match1=match1+1; end
        if dec_g2_1==dec_g2_2; match2=match2+1; end
    end
    runs(run,:)=[match1/num_trials,match2/num_trials];
end
toc
boxplot(runs)