% cross-param PP L-NL
% assesses PP of L model to NL model for different (gamma,h) pairs.

clear
parpool([12,80])
tic
rng('shuffle')
nsd=1; % noise
ntrials=100000;
gammas=0:0.1:10; num_gammas=length(gammas); 
hs=0:0.1:2.5; num_h=length(hs);

% ------------------Get a bank of clicks data-----------------------------%

db = '/scratch/adrian/validation2.h5';%'../data/validation2.h5';%'/home/adrian/validation1.h5';
[trials,envt]=get_trials(db,ntrials);

% remove problematic trials
bad_trials=[93737,207048,229626,272270,555142,631387,666886,774387,811388,961053];
to_remove=bad_trials(bad_trials<=ntrials); num_bad=length(to_remove);
if num_bad
    trials(:,to_remove)=[];
    ntrials=ntrials-length(to_remove);
end

left_clicks=trials(1,:);
right_clicks=trials(2,:);

high_rate=20; low_rate=5; k=log(high_rate/low_rate);

PP=zeros(num_h,num_gammas);
for idx_h=1:num_h
    for idx_g=1:num_gammas
        
        h=hs(idx_h); g=gammas(idx_g);
        
        match_count=0;
        
        parfor trn=1:ntrials
            lst=left_clicks{trn}; rst=right_clicks{trn};
            %         [lst, rst]=trials{1:2,trn};
            
            total_clicks = length(lst)+length(rst);
            
            % generate decisions with stoch nonlinear model
            dec_h = decide_AR(2,lst,rst,NaN,h,0,NaN,...
                normrnd(k, nsd, [total_clicks, 1]));
            
            % generate decisions with stoch linear model
            dec_g = gauss_noise_lin_decide(lst, rst, g, k, nsd, 0);
            
            % flip a coin if any decision was 0
            if dec_h == 0; dec_h = sign(rand-0.05); end
            if dec_g == 0; dec_g = sign(rand-0.05); end
            
            if dec_h==dec_g
                match_count=match_count+1;
            end
        end
        PP(idx_h,idx_g)=match_count/ntrials;
    end
end
toc
savefile=['/home/adrian/joint_PP_ntrials_',num2str(ntrials),...
    '_noise_',num2str(nsd),'.mat'];
save(savefile,'PP','gammas','hs','ntrials','-v7.3')
