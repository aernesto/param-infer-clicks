function trials=fetch_trials(filename,trial_range)
% ARGS:
%   filename        full path to .h5 file
%   trial_range     1-by-2 vector describing a range of integers (endpoints
%                   included)
% RETURNS:
%   trials          3-by-N cell, where N=trial_range(2)-trial_range(1)+1. 
%                   First two rows of the cell contain click times for left
%                   and right streams. Third row contains change point 
%                   times. So, to access left and right click streams from
%                   trial M, as col vectors, type: 
%                   [left,right]=trials{1:2,M};
file_info = h5info(filename);
group_name = file_info.Groups.Name;
   
trials = h5read(filename, [group_name,'/trials']);

trials_dset_name=[group_name,'/trials'];
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset
if tot_num_trials < trial_range(2) || trial_range(1)<1
    err('trial_range out of bounds')
else
    trials=trials(:,trial_range(1):trial_range(2));
end
end