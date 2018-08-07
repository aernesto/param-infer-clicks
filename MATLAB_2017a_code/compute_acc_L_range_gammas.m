% compute acc of L model for a range of gamma values on validation2.h5
clear
tic
rng('shuffle')
parpool([12,80])
ntrials=1000000;%10000;takes 24 sec on laptop
gammas=0:.05:20; num_gammas=length(gammas);
db='/scratch/adrian/validation2.h5';%'../data/validation2.h5';

%1. -----------------------get trials-------------------------------------%

[trials,envt]=get_trials(db,ntrials);
high_rate=20; low_rate=5; k=log(high_rate/low_rate);

%2. --------------------------compute acc --------------------------------%

Correct= zeros(ntrials,num_gammas);
nsd = 1; % noise
parfor trn=1:ntrials
    [lst,rst]=trials{:,trn};
    
    % generate decisions with stoch linear model
    if isempty(rst)
        rst = -Inf; right_noise=0; 
    else
        right_noise=normrnd(k,nsd,...
            [length(rst),num_gammas]);
    end
    if isempty(lst)
        lst = -Inf; left_noise=0; 
    else
        left_noise=normrnd(k,nsd,...
            [length(lst),num_gammas]);
    end
    
    dec = sign(sum(right_noise.*exp(rst*gammas),1)...
                    -sum(left_noise.*exp(lst*gammas),1));
    num_zeros = length(find(~dec));
    dec(dec==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = dec == envt(2,trn);
end
acc = sum(Correct,1)/ntrials;

%3. -----------------store and save accuracy------------------------------% 
%plot(gammas,acc)
savefile='/home/adrian/acc_L_range_gammas.mat';
save(savefile,'acc','gammas')
toc