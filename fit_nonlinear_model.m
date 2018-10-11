function [posterior, point_estimate]=fit_nonlinear_model(ref_model, dbname,...
    disc_prior,trial_range,ndiscount,npart,shuffle_db)
% fits the stochastic linear model to data
% ARGUMENTS:
%   ref_model   -- string. either 'lin' or 'nonlin'
%   dbname      -- string. full path to .h5 db file
%   disc_prior  -- 1-by-2 vector describing the endpoints of the support
%                   for the prior over the discounting parameter
%   trial_range -- 1-by-2; trials to use for the fit 
%   ndiscount   -- number of discounting param values to use for likelihood
%   npart       -- number of particles to use to estimate likelihood
%   shuffle_db  -- boolean. If true, trials are randomly permuted
% WARNING: The shuffle_db=true flag doesn't currently work
% NOTES:
%   Calls: llh2density() & lhd_nonlin_sing_tr_gauss_clicks()

% sample space (vector of samples to use. 1 sample = 1 discounting rate)
hstart=disc_prior(1);hend=disc_prior(2);
hs=linspace(hstart, hend, ndiscount)';
dh=(hend-hstart)/(ndiscount-1);  % step between consecutive samples

% database info (where the clicks data and other parameter values reside)
params=fetch_params(dbname, ref_model);
trials = fetch_trials(dbname,trial_range); % clicks data

% extract reference decisions into row vector:
ref_decisions = fetch_model_responses(dbname,trial_range,ref_model);

tic

if shuffle_db
    % shuffle trial order
    err('shuffle feature not enabled yet')
%    rng('shuffle')
%    trials = trials(:,randperm(tot_trials_db));
end

llh = zeros(ndiscount,1);
num_trials=trial_range(2)-trial_range(1)+1;
parfor trn=1:num_trials
    [lst, rst]=trials{1:2,trn};
    total_clicks = length(lst)+length(rst);
    
    % read synthetic decision from db
    synthetic_decision = ref_decisions(trn);
    
    % for random numbers generations to come
    rng('shuffle')
    
    % flip a coin if decision was 0
    if synthetic_decision == 0
        synthetic_decision = sign(rand-0.05);
    end
    
    % generate upfront all the Gaussian r.v. needed
    Gaussian_bank = normrnd(params.kappa, params.noise,...
        [npart, total_clicks, ndiscount]);
    
    % compute log-lh of each sample, for the 2 model pairs
    lhv=lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart, lst,...
        rst, params.T, hs, 0, Gaussian_bank);
    lhv(lhv<eps) = eps;
    lhd=log(lhv);
    llh=llh+lhd;
end

% shift log-likelihood up to avoid numerical errors
[max_val,idx1]=max(llh);
point_estimate=hs(idx1);

llh=llh+abs(max_val);

posterior=llh2density(llh,dh);           % convert log-lh to density

% plots for debugging
% subplot(2,1,1); plot(hs,llh)
% subplot(2,1,2); plot(hs,posterior)
fprintf('elapsed time to fit nonlinear model to %s model\n',ref_model)
toc
end
