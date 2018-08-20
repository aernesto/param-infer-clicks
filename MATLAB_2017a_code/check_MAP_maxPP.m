% attempt to estimate whether our MAP estimate is equivalent to maximizing
% PP
clear
rng('shuffle')

% sample space (vector of samples to use. 1 sample = 1 discounting rate)
gstart=0;gend=40; ndiscount=200;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);  % step between consecutive samples
ntrials=100000;      % number of trials to use in fitting procedure (for each block)

% database info (where the clicks data and other parameter values reside)
filename = '../data/validation2.h5';%'../data/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');   % hazard rate
true_g = 10; % reference discounting rate for linear model 
T = h5readatt(filename, info_dset_name,'T');        % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); % click rates
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
k=log(high_rate/low_rate);                          % jump in evidence at clicks
all_trials = h5read(filename, [group_name,'/trials']); % clicks data
tot_trials_db = size(all_trials,2);                 % total number of trials in DB

all_trials = all_trials(1:2,:);

nsd=2; % Gaussian noise applied to click height

nruns=1; % number of blocks of trials. MSE is computed across blocks
mses=0;  % MSE computed as running average
infs=zeros(1,ntrials);  % will store where log-likelihood = -Inf

tic

llh = zeros(ndiscount,1);
display_period=1000;
maps=zeros(1,floor(ntrials/display_period)); maps_idx=1;
trns=maps;
for trn=1:ntrials
    [lst, rst]=all_trials{:,trn};
    total_clicks = length(lst)+length(rst);
    refdec_noise = normrnd(k,nsd, [total_clicks,1]);
    
    % generate synthetic decision with nonlinear model
    synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);
    
    % generate synthetic decision with linear model
    %synthetic_decision = gauss_noise_lin_decide(lst, rst, true_g, k, nsd, 0);
    
    % flip a coin if decision was 0
    if synthetic_decision == 0
        synthetic_decision = sign(rand-0.05);
    end
    
    % compute log-lh of each sample
    % PB HERE!
    lhd=lhd_lin_sing_tr_gauss_clicks(synthetic_decision,...
        nsd, k, T, lst', rst', gs'); % already log-likelihood
    infs(trn)=sum(lhd == -Inf);      % check whether log-lh is -Inf
    llh=llh+lhd;
    if ~mod(trn,display_period)
        [~,map_idx]=max(llh);
        maps(maps_idx)=gs(map_idx);
        trns(maps_idx)=trn;
        maps_idx=maps_idx+1;
        density=llh2density_AR(llh,dg);                % convert log-lh to density
        mses=mses+dg*sum(((gs-true_g).^2).*density');  % running average
        %fprintf('trial %d \n',trn)
        %fprintf('mse %f \n',mses);
    end
end


figure()
plot(infs)
xlabel('trial number in block')
ylabel('number of -Inf in log-probs vector')

figure()
plot(trns, maps,'*','MarkerSize',14)
xlabel('trial number')
ylabel('MAP')
title('L-L noise=2 true \gamma=10')

mses=mses/nruns;

toc

%fname=['mses',num2str(nruns),'runs',num2str(ntrials),'trials'];
%save(['/home/adrian/tosubmit_home/',fname,'.mat'],'mses')