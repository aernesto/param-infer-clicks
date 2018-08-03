% estimate best gamma and best h for stochastic models
% only section 1 and then the appropriate section should be run

%% 1. get the trials
clear
filename = '../data/validation2.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
% nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
%lin_decision_dset_name=[group_name,'/decision_lin'];
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
%tot_num_trials = 500;
% h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate'); 
k=log(high_rate/low_rate); % kappa = jump size
%true_g = h5readatt(filename, lin_decision_dset_name, 'best_gamma');
all_trials = h5read(filename, [group_name,'/trials']);

% to get the start and end states of the environment for each trial, we do:
all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);
% row 1 of the cell corresponds to starting state of the environment
% row 2 of the cell corresponds to ending state of the environment
% a value of 1 indicates H+, -1 indicates H-

% example of extracting click times for left and right streams (col vecs)
% for trial 50
% [left_train, right_train] = all_trials{1:2,50};  

%% 2. run deterministic linear model on trials 

tic
gammas=6.8:0.005:7.1; num_gammas=length(gammas);
Correct= zeros(tot_num_trials,num_gammas);
for trn=1:tot_num_trials
    [lst, rst]=all_trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
    
    % generate decisions with deterministic linear model
    if isempty(rst); rst = -Inf; end
    if isempty(lst); lst = -Inf; end
    
    deter_dec = sign(sum(exp(rst*gammas),1)-sum(exp(lst*gammas),1)); 
    num_zeros = length(find(~deter_dec));
    deter_dec(deter_dec==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = deter_dec == all_envt(2,trn);
end
Acc = sum(Correct,1)/tot_num_trials;
[M,I]=max(Acc);
gammas(I)
plot(gammas,Acc)
toc

%% 3. run stochastic linear model on trials
rng('shuffle')
tic
gammas=5:0.1:7; num_gammas=length(gammas);
Correct= zeros(tot_num_trials,num_gammas);
for trn=1:tot_num_trials
    [lst, rst]=all_trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
    
    nsd = 2; % noise
    % generate decisions with stoch linear model
    if isempty(rst)
        rst = -Inf; right_noise=0; 
    else
        right_noise=normrnd(log(high_rate/low_rate),nsd,...
            [length(rst),num_gammas]);
    end
    if isempty(lst)
        lst = -Inf; left_noise=0; 
    else
        left_noise=normrnd(log(high_rate/low_rate),nsd,...
            [length(lst),num_gammas]);
    end
    
    deter_dec = sign(sum(right_noise.*exp(rst*gammas),1)...
                    -sum(left_noise.*exp(lst*gammas),1));
    num_zeros = length(find(~deter_dec));
    deter_dec(deter_dec==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = deter_dec == all_envt(2,trn);
end
Acc = sum(Correct,1)/tot_num_trials;
[M,I]=max(Acc);
gammas(I)
plot(gammas,Acc)
toc

%% 4. run stochastic nonlinear model on trials
rng('shuffle')
tic
h=0.3:0.01:.7; num_h=length(h);
Correct= zeros(tot_num_trials,num_h);
tot_num_trials = 100000;
nsd = 1; % noise

for trn=1:tot_num_trials
    [lst, rst]=all_trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
        
    % generate decisions with stoch nonlinear model
    decision_nonlin = decide_AR(T,...
        lst, rst, NaN, h', 0, NaN, normrnd(k, nsd, [total_clicks, 1]))';
   
    % flip a coin if any decision was 0
    num_zeros = length(find(~decision_nonlin));
    decision_nonlin(decision_nonlin==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = decision_nonlin == all_envt(2,trn);
end
Acc = sum(Correct,1)/tot_num_trials;
[M,I]=max(Acc);
h(I)
plot(h,Acc)
toc