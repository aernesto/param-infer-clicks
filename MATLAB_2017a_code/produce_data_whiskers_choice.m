%1. Get trials data
%2. Compute decisions from linear and stochastic true models
%3. loop over 500 blocks (fix block size)
    % compute and store decision of 4 combinations
%4. compute and plot % match

clear

tic


%1. --------------Get trials data-----------------------------------------%

filename = '../data/validation1.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
% nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
lin_decision_dset_name=[group_name,'/decision_lin'];
trial_info = h5info(filename, trials_dset_name);
%tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
tot_num_trials = 500;
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





%2. ----Compute decisions from linear and stochastic true models----------%

% left col for dec of lin, right col for dec of nonlin
decisions = zeros(tot_num_trials,2);

k=log(high_rate/low_rate);
nsd=1; % Gaussian noise applied to click height
rng(1) % for reference choice data to be reproducible
for trn=1:tot_num_trials
    [lst, rst]=all_trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
    refdec_noise = normrnd(k,nsd, [total_clicks,1]);
    
    % generate synthetic decision with nonlinear model
    synthetic_decision_nonlin = decide_AR(T, lst, rst, k, 1, 0, nsd,...
        refdec_noise);
    
    % generate synthetic decision with linear model
    synthetic_decision_lin = gauss_noise_lin_decide(lst, rst, true_g, k,...
        nsd, 0);
    
    % flip a coin if decision was 0
    if synthetic_decision_nonlin == 0
        synthetic_decision_nonlin = sign(rand-0.05);
    end
    if synthetic_decision_lin == 0
        synthetic_decision_lin = sign(rand-0.05);
    end
    
    % store decisions
    decisions(trn,:)=[synthetic_decision_lin, synthetic_decision_nonlin];
end
rng('shuffle')




%3&4. --------------------------------------------------------------------%
%loop over 500 blocks (fix block size)
    % compute and store decision of 4 combinations

nblocks = 500;  % number of blocks = number of points for single whisker
block_size=500; % block size = number of trials used for the fits
match = zeros(nblocks,4); % percentage of decision match
                    % COL 1 = LL
                    % COL 2 = L-NL
                    % COL 3 = NL-NL
                    % COL 4 = NL-L
                    
                    
                    
% load discounting parameter values for fitted models
% after the two load commands below, the variables modes_<mfmd> are in the
% workspace, where mf and md are one of lin and nonlin.
load(['../data/mse_nonlin_fig4_iteration2_',num2str(block_size),...
    'trials.mat'])
load(['../data/mse_lin_fig4_iteration2_',num2str(block_size),'trials.mat'])


for block=1:nblocks
    fit_decisions=zeros(tot_num_trials,4);  % col order is same as match 
        
    % Naming convention below
    % nonlin_1 means NL-NL
    % nonlin_2 means NL-L
    % lin_1    means L-L
    % lin_2    means L-NL
    h_nonlin_1 = modes_nonlinnonlin(block);
    h_nonlin_2 = modes_nonlinlin(block);
    g_lin_1 = modes_linlin(block);
    g_lin_2 = modes_linnonlin(block);
    
    for trn=1:tot_num_trials
        [lst, rst]=all_trials{1:2,trn};
        total_clicks = length(lst)+length(rst);
        nonlin_noise_1 = normrnd(k,nsd, [total_clicks,1]);
        nonlin_noise_2 = normrnd(k,nsd, [total_clicks,1]);
        % generate decisions with nonlinear models
        dec_nonlin_1 = decide_AR(T, lst, rst, k, h_nonlin_1, 0, nsd,...
            nonlin_noise_1);
        dec_nonlin_2 = decide_AR(T, lst, rst, k, h_nonlin_2, 0, nsd,...
            nonlin_noise_2);
        
        % generate decisions with linear models
        dec_lin_1 = gauss_noise_lin_decide(lst, rst, g_lin_1, k,...
            nsd, 0);
        dec_lin_2 = gauss_noise_lin_decide(lst, rst, g_lin_2, k,...
            nsd, 0);
        
        % flip a coin if decision was 0
        if dec_nonlin_1 == 0
            dec_nonlin_1 = sign(rand-0.05);
        end
        if dec_lin_1 == 0
            dec_lin_1 = sign(rand-0.05);
        end
        if dec_nonlin_2 == 0
            dec_nonlin_2 = sign(rand-0.05);
        end
        if dec_lin_2 == 0
            dec_lin_2 = sign(rand-0.05);
        end
        
        % store decisions
        fit_decisions(trn,:)=...
            [dec_lin_1,dec_lin_2,dec_nonlin_1,dec_nonlin_2];
    end
    for pair = 1:4
        if ismember(pair,[1,4])
            ref = 1;
        else
            ref = 2;
        end
        match(block, pair) = sum(fit_decisions(:,pair)==...
            decisions(:,ref)); 
    end
end
match = match / tot_num_trials;
toc
save(['../data/choice_match_',num2str(block_size),'_2.mat'])