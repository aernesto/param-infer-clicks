% performs the histogram test: https://paper.dropbox.com/doc/Battery-of-tests-for-model-fits-UkWXwtrYyULoGlxOlPTE8#:h2=Find-out-why-this-code-produce
clear
rng('shuffle')
lw=3;   % line width for plots
fs=25;fs2=20;  %font size for plots
ndiscount=2; % number of discounting parameter values to try
% sample space (vector of samples to use. 1 sample = 1 discounting rate)
gstart=0.5;gend=25;
gs=linspace(gstart, gend, ndiscount);
dg=(gend-gstart)/(ndiscount-1);  % step between consecutive samples
ntrial_vec=1;      % number of trials to use in fitting procedure 
                             %(for each block)
npart=800; % number of particles

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

tic
ntrials=[33,99];
for iii=1:2
    trn=ntrials(iii);
    [lst, rst]=all_trials{:,trn};
    total_clicks = length(lst)+length(rst);
    refdec_noise = normrnd(k,nsd, [total_clicks,1]);
    
    % generate synthetic decision with nonlinear model
    %synthetic_decision = decide_AR(T, lst, rst, k, true_h, 0, nsd, refdec_noise);
    
    % generate synthetic decision with linear model
    synthetic_decision_lin = gauss_noise_lin_decide(lst, rst, true_g, k, nsd, 0);
    
    % flip a coin if decision was 0
    %if synthetic_decision == 0
     %   synthetic_decision = sign(rand-0.05);
    %end
    if synthetic_decision_lin == 0
        synthetic_decision_lin = sign(rand-0.05);
    end
    
    % get sufficient statistics 
    linlin_theo=suff_stats_lhd_lin(synthetic_decision_lin,...
        nsd, k, T, lst', rst', gs');
    %linnonlin_theo=suff_stats_lhd_lin(synthetic_decision,...
     %   nsd, k, T, lst', rst', gs'); 
    
    % run particle filter
    % generate upfront all the Gaussian r.v. needed
    Gaussian_bank = normrnd(k, nsd, [npart, total_clicks, ndiscount]);
    landing_heights = hist_part_filt(npart,...
                                     lst,...
                                     rst,...
                                     T,...
                                     gs,...
                                     0,...
                                     Gaussian_bank);
    figure()
    for ii=1:2
        ax=subplot(2,1,ii);
        histogram(landing_heights(:,ii),'Normalization','pdf')
        xrange=xlim();x=xrange(1):.01:xrange(2);
        yrange=ylim;
        hold on
        plot(x,normpdf(x,linlin_theo(1,ii),linlin_theo(2,ii)),...
            'LineWidth',4)
        hold off
        xlim(xrange)
        ylim(yrange)
        ax.FontSize=20;
    end
end

toc