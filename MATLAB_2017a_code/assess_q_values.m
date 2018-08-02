% assess the distribution of the Bernoulli parameter q across clicks trials
% for both linear and nonlinear model
% Parameter q is defined as follows: it is the probability that the
% decision maker answers 1, given the clicks stimulus
% I assume all clicks stimuli have equal probability and only consider the
% histogram of q-values
clear
tic
%1. -------------Get a bank of clicks data---------------------------------

filename = '../data/validation1.h5';%'../data/S3lr5h1T2tr10000sp1000.h5';%'../data/validation1.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
% nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
lin_decision_dset_name=[group_name,'/decision_lin'];
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
%tot_num_trials = 1000;
npart=800;
% h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate'); 
true_g = h5readatt(filename, lin_decision_dset_name, 'best_gamma');
all_trials = h5read(filename, [group_name,'/trials']);
% all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);

% example of extracting click times for left and right streams (col vecs)
% for trial 50
% [left_train, right_train] = all_trials{1:2,50};  


%2. -------------For each trial, compute the q-value-----------------------

k=log(high_rate/low_rate);
nsd=1; % Gaussian noise applied to click height
q=zeros(1,tot_num_trials); % store q-values

for trn=1:tot_num_trials
    [lst, rst]=all_trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
    
    % for linear model, uncomment the following
    %q(trn)=exp(lhd_lin_sing_tr_gauss_clicks(1,nsd,k,T,lst',rst',true_g));
    
    
    noise_bank = normrnd(k, nsd, [npart, total_clicks, 1]);
    % for nonlinear model, uncomment the following
    q(trn)=lhd_AR(1, npart, lst, rst, 2, k, 1, 0, nsd, noise_bank);
end

nbins=100;  % actual number of bins for histogram and density
qbins=linspace(0,1,nbins+1);  % really bin edges as needed by histogram()
fig=figure();
fig.Visible='off';
hist=histogram(q,'Normalization','pdf','NumBins',nbins,'BinEdges',qbins);
save('../data/nonlin_q_density_10000trials.mat','hist', 'npart','tot_num_trials','nsd')
toc