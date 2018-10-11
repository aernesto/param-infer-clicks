function correct=fetch_correct_responses(filename,trial_range)
% ARGS:
%   filename        full path to .h5 file
%   trial_range     1-by-2 vector describing a range of integers (endpoints
%                   included)
% RETURNS:
%   correct         1-by-N vector, where N=trial_range(2)-trial_range(1)+1. 
%                   Vector contains -1 at entry n whenever correct response
%                   on trial n was -1 (i.e. H-). +1 codes for H+.
file_info = h5info(filename);
group_name = file_info.Groups.Name;
trials_dset_name=[group_name,'/trials'];
trial_info = h5info(filename, trials_dset_name);
all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);

tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
if tot_num_trials < trial_range(2) || trial_range(1)<1
    err('trial_range out of bounds')
else
    correct=all_envt(2,trial_range(1):trial_range(2));
end
end