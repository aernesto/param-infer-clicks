% compute acc of NL model for a range of h values on validation2.h5
clear
tic
rng('shuffle')
parpool([12,80])
ntrials=1000000;%10000;takes 24 sec on laptop
hs=0:.01:5; num_h=length(hs);
db='/scratch/adrian/validation2.h5';%'../data/validation2.h5';

%1. -----------------------get trials-------------------------------------%

[trials,envt]=get_trials(db,ntrials);
high_rate=20; low_rate=5; k=log(high_rate/low_rate);

%2. -------------------------compute acc----------------------------------%

Correct= zeros(ntrials,num_h);
nsd = 1; % noise
parfor trn=1:ntrials
    [lst,rst]=trials{:,trn};
    total_clicks = length(lst)+length(rst);
        
    % generate decisions with stoch nonlinear model
    dec = decide_AR(2,...
        lst, rst, NaN, hs', 0, NaN, normrnd(k, nsd, [total_clicks, 1]))';
   
    % flip a coin if any decision was 0
    num_zeros = length(find(~dec));
    dec(dec==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = dec == envt(2,trn);
end
acc = sum(Correct,1)/ntrials;

%3. -----------------store and save accuracy------------------------------% 
%plot(hs,acc)
savefile='/home/adrian/acc_NL_range_h.mat';
save(savefile,'acc','hs')
toc