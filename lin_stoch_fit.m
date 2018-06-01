clear
rng('shuffle')
lw=3; % line width for plots
fs=15;%font size for plots
ndiscount=800; % number of discounting parameter values to try
%hstart=0;hend=10; % range should be large enough for normalization
%hs=linspace(hstart,hend,ndiscount)'; % values of h to try
%dh=(hend-hstart)/(ndiscount-1);
gstart=0;gend=40;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);
ntrials=1;
filename = 'data/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');  % hazard rate
true_g = h5readatt(filename, [group_name, '/decision_lin'], 'best_gamma');
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
k=log(high_rate/low_rate);
all_trials = h5read(filename, [group_name,'/trials']);
tot_trials_db = size(all_trials,2);
% shuffle trial order
all_trials = all_trials(1:2,randperm(tot_trials_db));
nsd=1.5; % Gaussian noise applied to click height

nruns=1;
mses=0;
infs=zeros(1,ntrials);
tic
for run=1:nruns
    % reshuffle trials at each run
    all_trials = all_trials(:,randperm(tot_trials_db));
    llh = zeros(ndiscount,1);
    parfor trn=1:ntrials
        [lst, rst]=all_trials{:,trn};
        total_clicks = length(lst)+length(rst);
        refdec_noise = normrnd(k,nsd, [1, total_clicks]);
        synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);
        % generate upfront all the Gaussian r.v. needed
        %Gaussian_bank = normrnd(k, nsd, [npart, total_clicks]);
        lhd=lhd_lin_sing_tr_gauss_clicks(synthetic_decision,...
            nsd, k, T, lst', rst', gs'); % already log-likelihood
        infs(trn)=sum(lhd == -Inf);
        llh=llh+lhd;
    end
    density=llh2density_AR(llh,dg);
    mses=mses+dg*sum(((gs-true_g).^2).*density');
end
sum(infs)
mses=mses/nruns
toc
fname=['mses',num2str(nruns),'runs',num2str(ntrials),'trials'];
%save(['/home/adrian/tosubmit_home/',fname,'.mat'],'mses')
plot(gs, density,'LineWidth',lw);
hold on
plot([true_g, true_g], [0,max(density)], 'r', 'LineWidth', lw)
hold off
ylabel('likelihood','FontSize',fs)
xlabel('gamma values','FontSize',fs)
title(['noise=',num2str(nsd)],'FontSize',fs)
