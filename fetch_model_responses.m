function responses=fetch_model_responses(filename,trial_range,model_type)
% ARGS:
%   filename        full path to .h5 file
%   trial_range     1-by-2 vector describing a range of integers (endpoints
%                   included)
%   model_type      either 'lin' or 'nonlin'. model to fetch responses from
% RETURNS:
%   responses       1-by-N vector, where N=trial_range(2)-trial_range(1)+1. 
%                   Vector contains -1 at entry n whenever response
%                   on trial n was -1 (i.e. H-). +1 codes for H+.
file_info = h5info(filename);
group_name = file_info.Groups.Name;

decision_dset_name=[group_name,'/decision_',model_type];

trials_dset_name=[group_name,'/trials'];
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset

if tot_num_trials < trial_range(2) || trial_range(1)<1
    err('trial_range out of bounds')
else
    % somehow specifying the start/count args of h5read doesn't work
    responses=h5read(filename, decision_dset_name);
    responses=responses(trial_range(1):trial_range(2));
end
end