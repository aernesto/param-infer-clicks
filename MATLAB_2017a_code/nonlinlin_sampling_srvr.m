clear
rng('shuffle')
lw=3; % line width for plots
ndiscount=200; % number of discounting parameter values to try
hstart=0;hend=10; % range should be large enough for normalization
hs=linspace(hstart,hend,ndiscount)'; % values of h to try
dh=(hend-hstart)/(ndiscount-1);
ntrials=50;
filename = '/home/adrian/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
true_gamma = h5readatt(filename, [group_name,'/decision_lin'],'best_gamma');
k=log(high_rate/low_rate);
all_trials = h5read(filename, [group_name,'/trials']);
tot_trials_db = size(all_trials,2);
% shuffle trial order
all_trials = all_trials(1:2,randperm(tot_trials_db));
npart = 800;
nsd=1.5; % Gaussian noise applied to click height

nruns=500;
mses=0;
tic
for run=1:nruns
    all_trials = all_trials(:,randperm(tot_trials_db));
    llh = zeros(ndiscount,1);
    parfor trn=1:ntrials
        [lst, rst]=all_trials{:,trn};
        total_clicks = length(lst)+length(rst);
        synthetic_decision = gauss_noise_lin_decide(lst, rst, true_gamma, k, nsd, 0);
        % generate upfront all the Gaussian r.v. needed
        Gaussian_bank = normrnd(k, nsd, [npart, total_clicks]);
        llh=llh+log(lhd_AR(synthetic_decision, npart, lst, rst, T, k, hs, 0, nsd, Gaussian_bank));
    end
    density=llh2density_AR(llh,dh);
    mses=mses+dh*sum(((hs-true_h).^2).*density);
end
mses=mses/nruns;
toc
fname=['mses',num2str(nruns),'runs',num2str(ntrials),'trials'];
save(['/home/adrian/tosubmit_home/',fname,'.mat'],'mses')
%plot(hs, exp(llh),'LineWidth',lw);
%ylabel('likelihood')
%xlabel('h values')
%title(['noise=',num2str(nsd)])
