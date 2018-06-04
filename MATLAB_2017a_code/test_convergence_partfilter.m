% checks convergence as fcn of nb of particles
clear
rng('shuffle')
lw=3; % line width for plots

% db 
dbname='data/S3lr5h1T2tr10000sp1000.h5'; % hdf5 filename
grp_name='/lr5hr20h1T2';
dsetname = [grp_name,'/trials'];
info_dset = [grp_name,'/trial_info'];
dec_dset = [grp_name,'/decision_nonlin']; % decision dataset

hr = double(h5readatt(dbname, info_dset, 'high_click_rate'));
lr = double(h5readatt(dbname, info_dset, 'low_click_rate'));
T = double(h5readatt(dbname, info_dset, 'T')); % trial duration
true_h = double(h5readatt(dbname, info_dset, 'h')); % true hazard rate
k=log(hr/lr); % kappa for mean jump size in LLR at click

hs=9;%linspace(0,5,50)'; % values of h to try
ntrials=1;
ncols=2; %nb of columns in db
trial_data = h5read(dbname, dsetname, [1 100], [ncols ntrials]);
lst = trial_data{1,1}; 
rst = trial_data{2,1}; 
parts = 100:200:2000;
noises = linspace(.5,1.5,3);
likelihoods = zeros(length(parts), length(noises));
tic
for kk = 1:length(noises)
    nsd=noises(kk);
    synthetic_decision = gauss_noise_nonlin_decide(T, lst, rst, k, true_h, 0, nsd);
    for jj = 1:length(parts)
        npart=parts(jj);            
        likelihoods(jj,kk)=lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart, lst, rst, T, k, hs, 0, nsd);
    end
    subplot(1,3,kk)
    plot(parts, likelihoods(:,kk),'LineWidth', 3)
    xlabel('nb of samples')
    ylabel('P(data|h)')
    title(['noise=',num2str(nsd)])
    ax=gca; ax.FontSize=16;
end
toc