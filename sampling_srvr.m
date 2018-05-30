clear
rng('shuffle')
lw=3; % line width for plots
ndiscount=50; % number of discounting parameter values to try
hs=linspace(0,5,ndiscount)'; % values of h to try
ntrials=50;
filename = '/home/adrian/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
k=log(high_rate/low_rate);
all_trials = h5read(filename, [group_name,'/trials']);
all_trials = all_trials(1:2,:);
npart = 500;
noises = linspace(.5,1.5,3);
nsd=noises(2);
llh = zeros(ndiscount,1);
tic
parfor trn=1:ntrials
    [lst, rst]=all_trials{:,trn};
    total_clicks = length(lst)+length(rst);
    refdec_noise = normrnd(k,nsd, [1, total_clicks]);
    synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);
    % generate upfront all the Gaussian r.v. needed
    Gaussian_bank = normrnd(k, nsd, [npart, total_clicks]);
    llh=llh+log(lhd_AR(synthetic_decision, npart, lst, rst, T, k, hs, 0, nsd, Gaussian_bank));
end
toc
save('/home/adrian/tosubmit_home/sampling1.mat','llh')
%plot(hs, exp(llh),'LineWidth',lw);
%ylabel('likelihood')
%xlabel('h values')
%title(['noise=',num2str(nsd)])
