% estimates posterior over h using particle filters
clear
lw=3; % line width for plots

% db 
dbname='data/S3lr2h1T2tr5Ksp10K.h5'; % hdf5 filename
grp_name='/lr2hr14h1T2';
dsetname = [grp_name,'/trials'];
info_dset = [grp_name,'/trial_info'];
dec_dset = [grp_name,'/decision_nonlin']; % decision dataset

hr = h5readatt(dbname, info_dset, 'high_click_rate');
lr = h5readatt(dbname, info_dset, 'low_click_rate');
T = h5readatt(dbname, info_dset, 'T'); % trial duration
true_h = h5readatt(dbname, info_dset, 'h'); % true hazard rate
k=log(hr/lr); % kappa for mean jump size in LLR at click

hs=linspace(0,10,10)'; % values of h to try
ntrials=1;
npart=100;
nsd=1;
ncols=2; %nb of columns in db
trial_data = h5read(dbname, dsetname, [1 1], [ncols ntrials]);
lst = trial_data{1,1}; % 2nd idx to change
rst = trial_data{2,1}; % 2nd idx to change
log_posterior=zeros(1,length(hs));
tic
for num_s = 1:length(hs)
    synthetic_decision = gauss_noise_nonlin_decide(T, lst, rst, k, true_h, 0, nsd);
    log_posterior(num_s)=log_posterior(num_s)+...
        log(lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart, lst, rst, T, k, hs(num_s), 0, nsd));
end
toc
plot(hs, log_posterior, 'LineWidth', lw)
%ylim([0,max(log_posterior)])
ylabel('unnormalized log-posterior')
xlabel('h fit')