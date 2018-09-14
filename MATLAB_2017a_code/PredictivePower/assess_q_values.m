% assess the distribution of the Bernoulli parameter q across clicks trials
% for both linear and nonlinear model
% Parameter q is defined as follows: it is the probability that the
% decision maker answers 1, given the clicks stimulus
% I assume all clicks stimuli have equal probability and only consider the
% histogram of q-values
clear
tic
%1. -------------Get a bank of clicks data---------------------------------

filename = '../data/validation2.h5';%'../data/S3lr5h1T2tr10000sp1000.h5';%'../data/validation1.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
% nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
lin_decision_dset_name=[group_name,'/decision_lin'];
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
%tot_num_trials = 1000;
% npart=800;
% h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate'); 
% true_g = h5readatt(filename, lin_decision_dset_name, 'best_gamma');
all_trials = h5read(filename, [group_name,'/trials']);
% all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);

% example of extracting click times for left and right streams (col vecs)
% for trial 50
% [left_train, right_train] = all_trials{1:2,50};  

left_clicks=all_trials(1,:);
right_clicks=all_trials(2,:);
%2. -------------For each trial, compute the q-value-----------------------

k=log(high_rate/low_rate);
nsd=1; % Gaussian noise applied to click height
gs=0:.2:1;num_gs=length(gs);
all_g=zeros(num_gs);
tot_pairs=num_gs*(num_gs+1)/2;
q1=zeros(tot_pairs,tot_num_trials); % store q-values
q2=q1;

parfor trn=1:tot_num_trials
    lst=left_clicks{trn}; rst=right_clicks{trn};

    % for linear model, uncomment the following
    for jj=1:tot_pairs
        q1(jj,trn)=exp(lhd_lin_sing_tr_gauss_clicks(1,nsd,k,T,lst',rst',...
            gs(mapgidx(jj,1))));
        q2(jj,trn)=exp(lhd_lin_sing_tr_gauss_clicks(1,nsd,k,T,lst',rst',...
            gs(mapgidx(jj,2))));
    end

end

% nbins=100;  % actual number of bins for histogram and density
% qbins=linspace(0,1,nbins+1);  % really bin edges as needed by histogram()
% fig=figure();
% fig.Visible='off';
% hist=histogram(q,'Normalization','pdf','NumBins',nbins,'BinEdges',qbins);
save(['../data/nonlin_q_density_noise_',num2str(nsd),'_10000trials.mat'],...
    'q1','q2','gs','tot_num_trials','nsd')
toc

function idx=mapgidx(linidx,qtype)
    
end
