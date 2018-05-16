% estimates posterior over h using particle filters
clear
lw=3; % line width for plots

% db 
dbname='/storage/adrian/data_S_2_5.h5'; % hdf5 filename
grp_name='/lr1hr6.46h1T2';
dsetname = [grp_name,'/trials'];
info_dset = [grp_name,'/trial_info'];
dec_dset = [grp_name,'/decision_nonlin']; % decision dataset

hr = double(h5readatt(dbname, info_dset, 'high_click_rate'));
lr = double(h5readatt(dbname, info_dset, 'low_click_rate'));
T = double(h5readatt(dbname, info_dset, 'T')); % trial duration
true_h = double(h5readatt(dbname, info_dset, 'h')); % true hazard rate
k=log(hr/lr); % kappa for mean jump size in LLR at click

hs=linspace(0,5,50)'; % values of h to try
ntrials=400;
npart=500;
nsd=1;
ncols=2; %nb of columns in db
trial_data = h5read(dbname, dsetname, [1 1], [ncols ntrials]);
log_posterior=zeros(1,length(hs));
tic
for num_trial = 1:ntrials
    lst = trial_data{1,num_trial}; 
    rst = trial_data{2,num_trial}; 
    for num_s = 1:length(hs)
        synthetic_decision = gauss_noise_nonlin_decide(T, lst, rst, k, true_h, 0, nsd);
        log_posterior(num_s)=log_posterior(num_s)+...
            log(lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart, lst, rst, T, k, hs(num_s), 0, nsd));
    end
end
toc
filename=['/home/adrian/tosubmit_home/logpost_',...
		'tr',...
		num2str(ntrials),...
		'sd',...
		num2str(nsd),...
		'part',...
		num2str(npart),...
		'sp',...
		num2str(length(hs)),...
		'.mat'];
save(filename,'log_posterior','-v7.3')
