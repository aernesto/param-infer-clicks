function acc=accuracy(model_type, model_params, dbname, trial_range)
% computes the choice accuracy of a given model on a given set of trials
% ARGS: 
%   model_type      Either 'lin' or 'nonlin'
%   params          struct with 2 fields:
%                       disc = Value of discounting parameter
%                       noise = STDEV of noise
%   dbname          full path to .h5 file
%   trial_range     must be an interval that fits within the db size
% RETURNS:
%   acc             accuracy value between 0 and 1

% fetch trials (throw error if trial_range out of bounds) 
trials=fetch_trials(dbname,trial_range);

% fetch correct responses (row vec)
correct_choices=fetch_correct_responses(dbname,trial_range);

% fetch necessary parameters
db_params=fetch_params(dbname,model_type);

% compute model's choices
num_trials=length(correct_choices);
model_choices=zeros(1,num_trials);

rng('shuffle') % change seed if reproducibility desired

% this loop is probably not efficient, but it will do for now
for tr=1:num_trials
    [left_train,right_train]=trials{1:2,tr};
    if strcmp(model_type,'lin')
        model_choices(tr)=decide_lin_noise(left_train, right_train,...
            model_params.disc, db_params.kappa, model_params.noise, 0);
    elseif strcmp(model_type,'nonlin')
        num_clicks=length(left_train)+length(right_train);
        noise_matrix=normrnd(db_params.kappa,model_params.noise,...
            [num_clicks,1]);
        model_choices(tr)=decide_nonlin_noise(db_params.T,left_train,...
            right_train, model_params.disc, 0, noise_matrix);
    else
        err('wrong model_type argument');
    end
end

% compute accuracy
acc=sum(model_choices==correct_choices)/num_trials;

end