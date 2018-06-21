% fits the stochastic linear model to data
clear
rng('shuffle')
lw=3;   % line width for plots
fs=25;fs2=20;  %font size for plots
ndiscount=800; % number of discounting parameter values to try
% sample space (vector of samples to use. 1 sample = 1 discounting rate)
gstart=0;gend=40;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);  % step between consecutive samples
ntrial_vec=100:100:500;      % number of trials to use in fitting procedure 
                             %(for each block)

% database info (where the clicks data and other parameter values reside)
filename = '../data/S3lr5h1T2tr10000sp1000.h5';%'/home/adrian/S3lr5h1T2tr10000sp1000.h5';
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

nsd=2; % Gaussian noise applied to click height

nruns=500; % number of blocks of trials. MSE is computed across blocks


tic
for ntrials=ntrial_vec
    mse_linlin=0; mse_linnonlin=0;  % MSE computed as running average
    % store the modes of posteriors
    modes_linlin=zeros(1,nruns); modes_linnonlin=modes_linlin;
    for run=1:nruns
        % shuffle trial order
        all_trials = all_trials(:,randperm(tot_trials_db));
        llh = zeros(ndiscount,2);%col1=linlin; col2=linnonlin
        parfor trn=1:ntrials
            [lst, rst]=all_trials{:,trn};
            total_clicks = length(lst)+length(rst);
            refdec_noise = normrnd(k,nsd, [total_clicks,1]);
            
            % generate synthetic decision with nonlinear model
            synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);
            
            % generate synthetic decision with linear model
            synthetic_decision_lin = gauss_noise_lin_decide(lst, rst, true_g, k, nsd, 0);
            
            % flip a coin if decision was 0
            if synthetic_decision == 0
                synthetic_decision = sign(rand-0.05);
            end
            if synthetic_decision_lin == 0
                synthetic_decision_lin = sign(rand-0.05);
            end
            
            % compute log-lh of each sample, for the 2 model pairs
            lhdl=lhd_lin_sing_tr_gauss_clicks(synthetic_decision_lin,...
                nsd, k, T, lst', rst', gs');
            lhdnl=lhd_lin_sing_tr_gauss_clicks(synthetic_decision,...
                nsd, k, T, lst', rst', gs'); % already log-likelihood

            llh=llh+[lhdl,lhdnl];
        end
        
        % shift log-likelihood up to avoid numerical errors
        [max_linlin,idx1]=max(llh(:,1));
        modes_linlin(run)=gs(idx1);
        [max_linnonlin,idx2]=max(llh(:,2));
        modes_linnonlin(run)=gs(idx2);
        llh(:,1)=llh(:,1)+abs(max_linlin);
        llh(:,2)=llh(:,2)+abs(max_linnonlin);

        density_lin=llh2density_AR(llh(:,1),dg);                % convert log-lh to density
        mse_linlin=mse_linlin+dg*sum(((gs-true_g).^2).*density_lin');  % running average
        density=llh2density_AR(llh(:,2),dg);                % convert log-lh to density
        mse_linnonlin=mse_linnonlin+dg*sum(((gs-true_g).^2).*density');  % running average
    end
    mse_linlin=(mse_linlin/nruns)/(true_g^2);
    mse_linnonlin=(mse_linnonlin/nruns)/(true_g^2);
    fprintf(',stoch,linlin,%d,%.10f\n,stoch,linnonlin,%d,%.10f\n',ntrials,...
        mse_linlin,ntrials,mse_linnonlin)
    fname=['err_lin_fig4_noise_',num2str(nsd),'_trials_',num2str(ntrials)];
    save(['../data/',fname,'.mat'],'mse_linnonlin',...
      'mse_linlin', 'modes_linlin', 'modes_linnonlin')
end

toc

