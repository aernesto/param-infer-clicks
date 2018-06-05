% fits the stochastic linear model to data
clear
rng('shuffle')
lw=3;   % line width for plots
fs=15;  %font size for plots
ndiscount=800; % number of discounting parameter values to try
% sample space (vector of samples to use. 1 sample = 1 discounting rate)
gstart=0;gend=40;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);  % step between consecutive samples
ntrials=1;      % number of trials to use in fitting procedure (for each block)

% database info (where the clicks data and other parameter values reside)
filename = '../data/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');   % hazard rate
true_g = h5readatt(filename,...                     % discounting rate for linear model 
    [group_name, '/decision_lin'], 'best_gamma');
T = h5readatt(filename, info_dset_name,'T');        % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); % click rates
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
k=log(high_rate/low_rate);                          % jump in evidence at clicks
all_trials = h5read(filename, [group_name,'/trials']); % clicks data
tot_trials_db = size(all_trials,2);                 % total number of trials in DB

all_trials = all_trials(1:2,:);

nsd=1.5; % Gaussian noise applied to click height

nruns=1; % number of blocks of trials. MSE is computed across blocks
mses=0;  % MSE computed as running average
%infs=zeros(1,ntrials);  % will store where log-likelihood = -Inf

tic
for run=1:nruns
    llh = zeros(ndiscount,1);
    for trn=1:ntrials
        [lst, rst]=all_trials{:,trn};
        total_clicks = length(lst)+length(rst);
        refdec_noise = normrnd(k,nsd, [total_clicks,1]);
        
        % generate synthetic decision with nonlinear model
        %synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);

        % generate synthetic decision with linear model
        synthetic_decision = gauss_noise_lin_decide(lst, rst, true_g, k, nsd, 0);
        
        % flip a coin if decision was 0
        if synthetic_decision == 0
            synthetic_decision = sign(rand-0.5);
        end

        % compute log-lh of each sample
        % PB HERE!
        lhd=lhd_lin_sing_tr_gauss_clicks(synthetic_decision,...
            nsd, k, T, lst', rst', gs'); % already log-likelihood
        %infs(trn)=sum(lhd == -Inf);      % check whether log-lh is -Inf
        llh=llh+lhd;
    end
    density=llh2density_AR(llh,dg);                % convert log-lh to density
    mses=mses+dg*sum(((gs-true_g).^2).*density');  % running average
end

%sum(infs)

mses=mses/nruns;

toc

%fname=['mses',num2str(nruns),'runs',num2str(ntrials),'trials'];
%save(['/home/adrian/tosubmit_home/',fname,'.mat'],'mses')
plot(gs, density,'LineWidth',lw);
hold on
plot([true_g, true_g], [0,max(density)], 'r', 'LineWidth', lw)
hold off
ylabel('likelihood','FontSize',fs)
xlabel('gamma values','FontSize',fs)
title(['noise=',num2str(nsd), 'synthet. dec=',num2str(synthetic_decision)],'FontSize',fs)
