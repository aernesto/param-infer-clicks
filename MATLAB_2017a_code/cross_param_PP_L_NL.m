% cross-param PP L-NL
% assesses PP of L model to NL model for different (gamma,h) pairs.

clear
tic
rng('shuffle')
nsd=1; % noise
ntrials=10;%100000
gammas=0:0.05:20; num_gammas=length(gammas); 
hs=0:0.01:5; num_h=length(hs);
num_pairs=num_gammas*num_h;

pairs=zeros(num_pairs,2);
for idx_h=1:num_h
    for idx_g=1:num_gammas
        num_pair=sub2ind([num_h,num_gammas], idx_h, idx_g);
        pairs(num_pair,:)=[gammas(idx_g),hs(idx_h)];
    end
end

% ------------------Get a bank of clicks data-----------------------------%

db = '../data/validation2.h5';%'/home/adrian/validation1.h5';
[trials,envt]=get_trials(db,ntrials);
high_rate=20; low_rate=5; k=log(high_rate/low_rate);


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
