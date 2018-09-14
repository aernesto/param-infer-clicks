% cross-param PP: three model pairs are possible. L-L; L-NL; NL-NL
% assesses PP for different (theta1,theta2) pairs.

clear
model_pair={'L','NL'};
parpool([12,80])
tic
rng('shuffle')
nsd=2; % noise
ntrials=1000000;

if ismember('L',model_pair)
    gammas=0:0.1:10; num_gammas=length(gammas); 
end
if ismember('NL',model_pair)
    hs=0:0.1:2.5; num_h=length(hs);
end

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

if strcmp('L',model_pair{1})
    thetas_1=gammas;
    num_theta_1=num_gammas;
elseif strcmp('NL',model_pair{1})
    thetas_1=hs;
    num_theta_1=num_h;
end
if strcmp('L',model_pair{2})
    thetas_2=gammas;
    num_theta_2=num_gammas;
elseif strcmp('NL',model_pair{2})
    thetas_2=hs;
    num_theta_2=num_h;
end

pair1=model_pair{1}; pair2=model_pair{2};

PP=zeros(num_theta_1,num_theta_2);
for idx_theta_1=1:num_theta_1
    
    if strcmp(pair1,pair2)
        theta_2_start=idx_theta_1; % only triangular matrix for LL and NLNL
    else
        theta_2_start=1;
    end
    
    for idx_theta_2=theta_2_start:num_theta_2
        
        theta_1=thetas_1(idx_theta_1); theta_2=thetas_2(idx_theta_2);
        
        match_count=0;
        
        parfor trn=1:ntrials
            lst=left_clicks{trn}; rst=right_clicks{trn};
            %         [lst, rst]=trials{1:2,trn};
            
            if strcmp('NL',pair1)
                % generate decisions with stoch nonlinear model
                total_clicks = length(lst)+length(rst);
                dec_1 = decide_AR(2,lst,rst,NaN,theta_1,0,NaN,...
                    normrnd(k, nsd, [total_clicks, 1]));
            else
                dec_1=gauss_noise_lin_decide(lst, rst, theta_1,k,nsd,0);
            end
            if strcmp('L',pair2)
                % generate decisions with stoch linear model
                dec_2 = gauss_noise_lin_decide(lst, rst, theta_2, k, nsd, 0);
            else
                % generate decisions with stoch nonlinear model
                total_clicks = length(lst)+length(rst);
                dec_2 = decide_AR(2,lst,rst,NaN,theta_2,0,NaN,...
                    normrnd(k, nsd, [total_clicks, 1]));
            end
            
            % flip a coin if any decision was 0
            if dec_1 == 0; dec_1 = sign(rand-0.05); end
            if dec_2 == 0; dec_2 = sign(rand-0.05); end
            
            if dec_1==dec_2
                match_count=match_count+1;
            end
        end
        PP(idx_theta_1,idx_theta_2)=match_count/ntrials;
    end
end
toc
savefile=['/home/adrian/joint_PP_',model_pair{1},model_pair{2},...
    '_ntrials_',num2str(ntrials),'_noise_',num2str(nsd),'.mat'];
save(savefile,'PP','model_pair','thetas_1','thetas_2','ntrials','-v7.3')
