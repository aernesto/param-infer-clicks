% compute acc of L model for a range of gamma values on validation2.h5
clear

server_version = true;

tic

rng('shuffle')

if server_version
    parpool([12,80])
    ntrials=1000000;
    db='/scratch/adrian/validation2.h5';
else
    ntrials=10000;
    db='../data/validation2.h5';
end

gammas=0:.05:20; num_gammas=length(gammas);

%1. -----------------------get trials-------------------------------------%

[trials,envt]=get_trials(db,ntrials);
high_rate=20; low_rate=5; k=log(high_rate/low_rate);

%2. --------------------------compute acc --------------------------------%

Correct= zeros(ntrials,num_gammas);
nsd = 1; % noise
nan_trial_count=0;
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
    % debug
    num_dec_nan=sum(isnan(dec));
    if num_dec_nan
        nan_trial_count=nan_trial_count+1;
        fprintf('trial %d, dec_nan=%d\n',trn, num_dec_nan)
    else
        num_zeros = length(find(~dec));
        dec(dec==0)=randsample([-1,1],num_zeros,true);
        Correct(trn, :) = dec == envt(2,trn);
    end
end
fprintf('total of %d NaN trials\n',nan_trial_count)
acc = sum(Correct,1)/(ntrials-nan_trial_count);

%3. -----------------store and save accuracy------------------------------% 
%plot(gammas,acc)
if server_version
    savefile='/home/adrian/acc_L_range_gammas.mat';
    save(savefile,'acc','gammas')
else
    plot(gammas,acc)
end
toc