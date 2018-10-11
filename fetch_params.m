function params=fetch_params(dbname, model_type)
% fetches parameters associated with the given database and decision model
% ARGS:
%   dbname          full path to .h5 db file
%   model_type      one of 'lin' or 'nonlin'. Specifies the model which
%                   produced the choice data present in the database
% RETURNS:
%   params          struct with parameter values. Field names should be
%                   self explanatory
filename = dbname;

%- Data structure 
% The data inside the .h5 file is stored inside a 'group'
% To find the group name, you can do the following
file_info = h5info(filename);
group_name = file_info.Groups.Name;

% The group itself contains 4 datasets: 'trials', 'trial_info',
% 'decision_lin' and 'decision_nonlin':
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
lin_decision_dset_name=[group_name,'/decision_lin'];

%- Read off parameters
% The dataset named '/trials' contains click times and environment state
% values for all the trials in the group.
% To find out how many trials are in the file, we do:
trial_info = h5info(filename, trials_dset_name);
params.tot_db_trials = trial_info.Dataspace.Size(2);  

% To get the values of more parameters regarding the data, we access the 
% attributes of the other datasets, as follows:
params.stim_hazard = h5readatt(filename, info_dset_name,'h');  
params.T = h5readatt(filename, info_dset_name,'T'); 
params.low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
params.high_rate = h5readatt(filename, info_dset_name,'high_click_rate'); 
params.kappa = log(params.high_rate/params.low_rate);
params.S = h5readatt(filename, info_dset_name,'S');  
params.model_type = model_type;
if strcmp(model_type,'lin')
    params.disc = h5readatt(filename, lin_decision_dset_name, 'disc_lin');
    params.noise= h5readatt(filename, lin_decision_dset_name, 'noise');
elseif strcmp(model_type,'nonlin')
    params.disc = h5readatt(filename, nonlin_decision_dset_name,...
        'disc_nonlin');
    params.noise= h5readatt(filename, nonlin_decision_dset_name, 'noise');
end

end