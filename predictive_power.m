function pp=predictive_power(target_model_type, target_model_params,...
    ref_model_type, dbname, trial_range)
% computes the predictive power (i.e. percent match in responses, of a 
% given model on a given set of trials
% ARGS: 
%   target_model_type      Either 'lin' or 'nonlin'
%   target_model_params    struct with 2 fields:
%                               disc = Value of discounting parameter
%                               noise = STDEV of noise
%   ref_model_type         either 'lin' or 'nonlin'
%   dbname          full path to .h5 file
%   trial_range     must be an interval that fits within the db size
% RETURNS:
%   acc             accuracy value between 0 and 1

% fetch trials (throw error if trial_range out of bounds) 
trials=fetch_trials(dbname,trial_range);

% fetch responses from reference model (row vec)
ref_model_choices=fetch_model_responses(dbname,trial_range,ref_model_type);

% fetch necessary parameters
db_params=fetch_params(dbname,target_model_type);

% compute target model's choices
num_trials=length(ref_model_choices);
model_choices=zeros(1,num_trials);

rng('shuffle') % change seed if reproducibility desired

% this loop is probably not efficient, but it will do for now
for tr=1:num_trials
    [left_train,right_train]=trials{1:2,tr};
    if strcmp(target_model_type,'lin')
        model_choices(tr)=decide_lin_noise(left_train, right_train,...
            target_model_params.disc, db_params.kappa, target_model_params.noise, 0);
    elseif strcmp(target_model_type,'nonlin')
        num_clicks=length(left_train)+length(right_train);
        noise_matrix=normrnd(db_params.kappa,target_model_params.noise,...
            [num_clicks,1]);
        model_choices(tr)=decide_nonlin_noise(db_params.T,left_train,...
            right_train, target_model_params.disc, 0, noise_matrix);
    else
        err('wrong model_type argument');
    end
end

% compute predictive power (i.e. percent match)
pp=sum(model_choices==ref_model_choices)/num_trials;

end