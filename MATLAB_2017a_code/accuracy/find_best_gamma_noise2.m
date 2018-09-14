% estimate best gamma and best h for stochastic models
clear
parpool([12,80])
rng('shuffle')
tic
tot_num_trials=1000000;
nsd=2;
gammas=3:.01:6; num_gammas=length(gammas);
h=0:.005:.05; num_h=length(h);

% ------------------Get a bank of clicks data-----------------------------%

filename = '/scratch/adrian/validation2.h5';%'../data/validation2.h5';%'/home/adrian/validation1.h5';
%filename = '../data/validation2.h5';
[all_trials,all_envt]=get_trials(filename,tot_num_trials);

% remove problematic trials
bad_trials=[93737,207048,229626,272270,555142,631387,666886,774387,811388,961053];
to_remove=bad_trials(bad_trials<=tot_num_trials); num_bad=length(to_remove);
if num_bad
    all_trials(:,to_remove)=[];
    tot_num_trials=tot_num_trials-length(to_remove);
end

left_clicks=all_trials(1,:);
right_clicks=all_trials(2,:);

high_rate=20; low_rate=5; k=log(high_rate/low_rate); T=2;

%------------------------ assess accuracy --------------------------------%

Correct= zeros(tot_num_trials,num_gammas);
Correct_h= zeros(tot_num_trials,num_h);

parfor trn=1:tot_num_trials
    lst=left_clicks{trn}; rst=right_clicks{trn};

    total_clicks = length(lst)+length(rst);
    
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
    
    lin_dec = sign(sum(right_noise.*exp(rst*gammas),1)...
                    -sum(left_noise.*exp(lst*gammas),1));
    num_zeros = length(find(~lin_dec));
    lin_dec(lin_dec==0)=randsample([-1,1],num_zeros,true);

    Correct(trn, :) = lin_dec == all_envt(2,trn);
    
   
    % generate decisions with stoch nonlinear model
    decision_nonlin = decide_AR(T,...
        lst, rst, NaN, h', 0, NaN, normrnd(k, nsd, [total_clicks, 1]))';
   
    % flip a coin if any decision was 0
    num_zeros = length(find(~decision_nonlin));
    decision_nonlin(decision_nonlin==0)=randsample([-1,1],num_zeros,true);

    Correct_h(trn, :) = decision_nonlin == all_envt(2,trn);
    
end
Acc = sum(Correct,1)/tot_num_trials;
[M,I]=max(Acc);
gammas(I)
%subplot(1,2,1)
%plot(gammas,Acc)

Acc_h = sum(Correct_h,1)/tot_num_trials;
[M_h,I_h]=max(Acc_h);
h(I_h)
%subplot(1,2,2)
%plot(h,Acc_h)
toc
README='this dset was produced with script find_best_gamma.m from commit 1d56e1d + 1';
fsave=['/home/adrian/tosubmit_home/acc_noise_',num2str(nsd),'.mat'];
save(fsave,'h','Acc_h','gammas','Acc','tot_num_trials','nsd','README')
