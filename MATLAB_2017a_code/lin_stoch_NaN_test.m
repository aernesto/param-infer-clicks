% Find out why fits of the stochastic linear model produces NaN data
% ref: https://paper.dropbox.com/doc/Battery-of-tests-for-model-fits-UkWXwtrYyULoGlxOlPTE8#:h2=Find-out-why-this-code-produce
clear
rng('shuffle')
lw=3;   % line width for plots
fs=25;fs2=20;  %font size for plots
ndiscount=800; % number of discounting parameter values to try
% sample space (vector of samples to use. 1 sample = 1 discounting rate)
gstart=0;gend=40;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);  % step between consecutive samples
ntrial_vec=400;      % number of trials to use in fitting procedure 
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

nsd=1; % Gaussian noise applied to click height

nruns=10; % number of blocks of trials. MSE is computed across blocks
running_min_linlin=zeros(nruns,ntrial_vec);
running_max_linlin=running_min_linlin;
running_min_linnonlin=zeros(nruns,ntrial_vec);
running_max_linnonlin=running_min_linnonlin;
tic
for ntrials=ntrial_vec
    mse_linlin=0; mse_linnonlin=0;  % MSE computed as running average
    
    for run=1:nruns
        % shuffle trial order
        all_trials = all_trials(:,randperm(tot_trials_db));
        llh = zeros(ndiscount,2);%col1=linlin; col2=linnonlin
        
        %%% LOOP OVER TRIALS
        for trn=1:ntrials
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
            
            % adding unscaled log-likelihoods
            llh=llh+[lhdl,lhdnl];
            
            % computing running extreme vals for debugging
            running_min_linlin(run, trn)=min(llh(:,1));
            running_min_linnonlin(run,trn)=min(llh(:,2));
            running_max_linlin(run, trn)=max(llh(:,1));
            running_max_linnonlin(run,trn)=max(llh(:,2));
            
        end
                
        % fitting to linear model
        density_lin=llh2density_AR(llh(:,1),dg);% convert log-lh to density
        
        toadd_lin=dg*sum(((gs-true_g).^2).*density_lin');
        
        mse_linlin=mse_linlin+toadd_lin;    % running average across blocks
        
        % fitting to nonlinear model
        density=llh2density_AR(llh(:,2),dg);    % convert log-lh to density
        
        toadd_nonlin=dg*sum(((gs-true_g).^2).*density');
        
        mse_linnonlin=mse_linnonlin+toadd_nonlin;% running average across blocks
        
        % DEBUG
        if sum(isnan(density_lin))
            fprintf('run %d, density_lin has %d NaNs\n',run, sum(isnan(density_lin)))
        end
        if sum(isnan(density))
            fprintf('run %d, density has %d NaNs\n',run, sum(isnan(density)))
        end
        if sum(isnan(toadd_lin))
            fprintf('run %d, toadd_lin has %d NaNs\n',run, sum(isnan(toadd_lin)))
        end
        if sum(isnan(toadd_nonlin))
            fprintf('run %d, toadd_nonlin has %d NaNs\n',run, sum(isnan(toadd_nonlin)))
        end
    end
    mse_linlin=(mse_linlin/nruns)/(true_g^2);
    mse_linnonlin=(mse_linnonlin/nruns)/(true_g^2);
    fprintf(',stoch,linlin,%d,%.10f\n,stoch,linnonlin,%d,%.10f\n',ntrials,...
        mse_linlin,ntrials,mse_linnonlin)
end


toc
figure()
subplot(1,2,1)
plot(1:ntrial_vec,running_min_linlin','LineWidth',lw)
title('running min linlin')
ax=gca;ax.FontSize=fs2;
subplot(1,2,2)
plot(1:ntrial_vec, running_min_linnonlin','LineWidth',lw)
title('running min linnonlin')
ax=gca;ax.FontSize=fs2;
figure()
subplot(1,2,1)
plot(1:ntrial_vec,running_max_linlin', 'LineWidth',lw)
title('running max linlin')
ax=gca;ax.FontSize=fs;
subplot(1,2,2)
plot(1:ntrial_vec, running_max_linnonlin','LineWidth',lw)
title('running max linnonlin')
ax=gca;ax.FontSize=fs;
% 
% fname=['mse_lin_fig4_iteration1_',num2str(ntrials),'trials'];
% save(['/home/adrian/tosubmit_home/',fname,'.mat'],'mse_linnonlin',...
%      'mse_linlin')