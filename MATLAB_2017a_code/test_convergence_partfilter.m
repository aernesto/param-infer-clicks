% checks convergence as fcn of nb of particles
clear
rng('shuffle')
lw=3; % line width for plots

% db 
dbname='../data/S3lr5h1T2tr10000sp1000.h5'; % hdf5 filename
grp_name='/lr5hr20h1T2';
dsetname = [grp_name,'/trials'];
info_dset = [grp_name,'/trial_info'];
dec_dset = [grp_name,'/decision_nonlin']; % decision dataset

hr = double(h5readatt(dbname, info_dset, 'high_click_rate'));
lr = double(h5readatt(dbname, info_dset, 'low_click_rate'));
T = double(h5readatt(dbname, info_dset, 'T')); % trial duration
true_h = double(h5readatt(dbname, info_dset, 'h')); % true hazard rate
k=log(hr/lr); % kappa for mean jump size in LLR at click

%%%% TO CHANGE
hs=.1;%linspace(0,5,50)'; % values of h to try
%%%%

ntrials=1;
ncols=2; %nb of columns in db

%%%%% TO CHANGE
trial_data = h5read(dbname, dsetname, [1 330], [ncols ntrials]);
%%%%%

lst = trial_data{1,1}; 
rst = trial_data{2,1}; 
nclicks=length(lst)+length(rst);
parts = 100:100:1600;
noises = 1;
likelihoods = zeros(length(parts), length(noises));
tic
for kk = 1:length(noises)
    nsd=noises(kk);
    
    % produce synthetic decision datum on trial stimulus
    %%%%% TO CHANGE
    noise1=normrnd(k, nsd, nclicks, 1);
    synthetic_decision = gauss_noise_nonlin_decide(T, lst, rst, k,...
        true_h, 0, nsd, noise1);
    %%%%%
    
    % run particle filter to estimate likelihood of decision for given hs
    for jj = 1:length(parts)
        npart=parts(jj);            
        %%%%% TO CHANGE
        noisebank=normrnd(k, nsd, npart, nclicks);
        likelihoods(jj,kk)=...
            lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart,...
            lst, rst, T, k, hs, 0, nsd, noisebank);
        %%%%%
    end
    
    % plot
    plot(parts, likelihoods(:,kk),'LineWidth', 3)
    hold on
    plot([800,800],[0,1],'r','LineWidth',2)
    hold off
    ylim([0,1])
    xlim([100,1600])
    ax=gca; 
    ax.YTick=[0,0.5,1];
    ax.XTick=[100,800,1600];
    ax.FontSize=16;
end
toc