% cross-param PP L
% assesses PP of L model to itself for different gamma pairs.
clear
rng('shuffle')
tic
%1. ---------------sanity check with same param---------------------------%
% Here, the goal is to see whether full numerical estimation of PP agrees
% with the semi-analytic way described here:
% https://paper.dropbox.com/doc/Stochastic-model-fitting--AJgGggz4kkvbqatj~ozkcoInAg-URuXW5PKABizdbMJB4Bkb#:uid=256420129419350336086512&h2=Predictive-power

gammas=0:.5:4; num_gammas=length(gammas);
nsd=1;
% ------------------Get a bank of clicks data-----------------------------%

db = '../data/validation2.h5';%'../data/validation1.h5';
file_info = h5info(db);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
% nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
%lin_decision_dset_name=[group_name,'/decision_lin'];
trial_info = h5info(db, trials_dset_name);
%tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
tot_num_trials = 20000;
% h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(db, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(db, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(db, info_dset_name,'high_click_rate'); 
k=log(high_rate/low_rate);
%true_g = 6.9;%6.7457;%h5readatt(filename, lin_decision_dset_name, 'best_gamma');
all_trials = h5read(db, [group_name,'/trials']);
% all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);

% example of extracting click times for left and right streams (col vecs)
% for trial 50
% [left_train, right_train] = all_trials{1:2,50};  


%----------------------------
% PP=zeros(size(gammas));
% for idx1=1:num_gammas
%     g=gammas(idx1);
%     match_count=0;
%     for trn=1:tot_num_trials
%         [lst, rst]=all_trials{1:2,trn};    
%         
%         dec_1 = gauss_noise_lin_decide(lst, rst, g, k, nsd, 0); 
%         if dec_1 == 0; dec_1 = sign(rand-0.05); end  % flip a coin if 0
%         
%         dec_2 = gauss_noise_lin_decide(lst, rst, g, k, nsd, 0);
%         if dec_2 == 0; dec_2 = sign(rand-0.05); end
%         
%         if dec_1==dec_2
%             match_count=match_count+1;
%         end
%     end
%     PP(idx1)=match_count/tot_num_trials;
% end
% toc
% plot(gammas,PP)

% ---------------------- end of sanity check -----------------------------%
PP=zeros(num_gammas);
for idx1=1:num_gammas
    for idx2=idx1:num_gammas
        g1=gammas(idx1); g2=gammas(idx2);
        match_count=0;
        for trn=1:tot_num_trials
            [lst, rst]=all_trials{1:2,trn};
            
            dec_1 = gauss_noise_lin_decide(lst, rst, g1, k, nsd, 0);
            if dec_1 == 0; dec_1 = sign(rand-0.05); end  % flip a coin if 0
            
            dec_2 = gauss_noise_lin_decide(lst, rst, g2, k, nsd, 0);
            if dec_2 == 0; dec_2 = sign(rand-0.05); end
            
            if dec_1==dec_2
                match_count=match_count+1;
            end
        end
        PP(idx1,idx2)=match_count/tot_num_trials;
    end  
end
toc
surf(PP)